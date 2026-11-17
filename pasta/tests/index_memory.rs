//! Reading corrupt BAI, CSI, tabix, and VCF indexes in a process that cannot allocate more than
//! 64 MiB at once.
//!
//! The counts in an index are untrusted, and every entry they count takes more memory than the
//! bytes it is read from. These indexes state counts that would need gigabytes, or hold more
//! empty reference sequences than fit under the limit. Reading them must fail with an error,
//! not abort the process on a failed allocation, whatever the machine's memory.

use std::{
    alloc::{GlobalAlloc, Layout, System},
    error::Error,
    io::{self, Write},
    ptr,
};

use pasta::{bam::bai, bgzf, csi, tabix, variant::index::Index};

/// The largest allocation this test binary makes.
const LIMIT: usize = 64 << 20;

/// The system allocator, failing every allocation larger than [`LIMIT`].
struct Limited;

// SAFETY: every call is passed on to the system allocator unchanged, or fails with null.
unsafe impl GlobalAlloc for Limited {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if layout.size() > LIMIT {
            ptr::null_mut()
        } else {
            // SAFETY: the caller upholds `GlobalAlloc::alloc`'s contract.
            unsafe { System.alloc(layout) }
        }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        if layout.size() > LIMIT {
            ptr::null_mut()
        } else {
            // SAFETY: the caller upholds `GlobalAlloc::alloc_zeroed`'s contract.
            unsafe { System.alloc_zeroed(layout) }
        }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: the caller upholds `GlobalAlloc::dealloc`'s contract.
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if new_size > LIMIT {
            ptr::null_mut()
        } else {
            // SAFETY: the caller upholds `GlobalAlloc::realloc`'s contract.
            unsafe { System.realloc(ptr, layout, new_size) }
        }
    }
}

#[global_allocator]
static ALLOCATOR: Limited = Limited;

/// The number of empty reference sequences in the indexes that hold too many: each takes at
/// most 8 bytes of index and far more than 32 bytes of memory, more than the limit in all.
const MANY: u32 = 1 << 21;

/// Compresses `data` in BGZF blocks.
fn bgzf_compress(data: &[u8]) -> Vec<u8> {
    let mut writer = bgzf::io::Writer::new(Vec::new());
    writer.write_all(data).expect("the data is compressed");
    writer.finish().expect("the data is compressed")
}

/// A CSI index (min shift 14, depth 5, no aux data) of `n_ref` reference sequences, followed by
/// `rest`.
fn csi(n_ref: u32, rest: &[u8]) -> Vec<u8> {
    let mut data = b"CSI\x01".to_vec();
    for n in [14, 5, 0, n_ref] {
        data.extend(n.to_le_bytes());
    }
    data.extend(rest);
    data
}

/// A tabix index of a VCF of `n_ref` reference sequences with no names, followed by `rest`.
fn tbi(n_ref: u32, rest: &[u8]) -> Vec<u8> {
    let mut data = b"TBI\x01".to_vec();
    // n_ref, format (VCF), col_seq, col_beg, col_end, meta ('#'), skip, l_nm.
    for n in [n_ref, 2, 1, 2, 0, u32::from(b'#'), 0, 0] {
        data.extend(n.to_le_bytes());
    }
    data.extend(rest);
    data
}

/// A BAI index of `n_ref` reference sequences, followed by `rest`.
fn bai(n_ref: u32, rest: &[u8]) -> Vec<u8> {
    let mut data = b"BAI\x01".to_vec();
    data.extend(n_ref.to_le_bytes());
    data.extend(rest);
    data
}

/// Returns whether `e` or an error it wraps is of the kind `kind`.
fn has_kind(e: &io::Error, kind: io::ErrorKind) -> bool {
    let mut source: Option<&(dyn Error + 'static)> = Some(e);
    while let Some(e) = source {
        if e.downcast_ref::<io::Error>()
            .is_some_and(|e| e.kind() == kind)
        {
            return true;
        }
        source = e.source();
    }
    false
}

/// Asserts that reading an index failed with an error of the kind `kind`.
fn assert_fails<T>(result: io::Result<T>, kind: io::ErrorKind) {
    match result {
        Ok(_) => panic!("the index was read"),
        Err(e) => assert!(has_kind(&e, kind), "{e:?}"),
    }
}

// Too many reference sequences to fit under the limit, each one empty: `n_bin` = 0, and in BAI
// and tabix `n_intv` = 0.

#[test]
fn bai_with_too_many_reference_sequences() {
    let data = bai(MANY, &vec![0; 8 * MANY as usize]);
    let result = bai::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::OutOfMemory);
}

#[test]
fn csi_with_too_many_reference_sequences() {
    let data = bgzf_compress(&csi(MANY, &vec![0; 4 * MANY as usize]));
    let result = csi::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::OutOfMemory);
}

#[test]
fn tabix_with_too_many_reference_sequences() {
    let data = bgzf_compress(&tbi(MANY, &vec![0; 8 * MANY as usize]));
    let result = tabix::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::OutOfMemory);
}

#[test]
fn vcf_csi_with_too_many_targets() {
    let result = Index::parse(&csi(MANY, &vec![0; 4 * MANY as usize]));
    assert_fails(result, io::ErrorKind::OutOfMemory);
}

#[test]
fn vcf_tbi_with_too_many_targets() {
    let result = Index::parse(&tbi(MANY, &vec![0; 8 * MANY as usize]));
    assert_fails(result, io::ErrorKind::OutOfMemory);
}

// A stated count with no entries after it: nothing is reserved for it, so the read fails where
// the data ends.

#[test]
fn bai_with_a_reference_sequence_count_past_the_data() {
    let data = bai(u32::MAX, &[]);
    let result = bai::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::UnexpectedEof);
}

#[test]
fn bai_with_a_bin_count_past_the_data() {
    let data = bai(1, &u32::MAX.to_le_bytes());
    let result = bai::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::UnexpectedEof);
}

#[test]
fn csi_with_a_reference_sequence_count_past_the_data() {
    let data = bgzf_compress(&csi(i32::MAX as u32, &[]));
    let result = csi::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::UnexpectedEof);
}

#[test]
fn csi_with_a_bin_count_past_the_data() {
    let data = bgzf_compress(&csi(1, &i32::MAX.to_le_bytes()));
    let result = csi::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::UnexpectedEof);
}

#[test]
fn tabix_with_a_reference_sequence_count_past_the_data() {
    let data = bgzf_compress(&tbi(i32::MAX as u32, &[]));
    let result = tabix::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::UnexpectedEof);
}

#[test]
fn tabix_with_a_bin_count_past_the_data() {
    let data = bgzf_compress(&tbi(1, &i32::MAX.to_le_bytes()));
    let result = tabix::io::Reader::new(&data[..]).read_index();
    assert_fails(result, io::ErrorKind::UnexpectedEof);
}
