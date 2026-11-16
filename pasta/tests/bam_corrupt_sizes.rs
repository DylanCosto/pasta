//! Reading BAM whose counts and sizes are corrupt.
//!
//! A count or size in a BAM header or record is untrusted: memory may grow only with the bytes
//! actually present. This binary's allocator refuses any single allocation larger than
//! [`MAX_ALLOCATION`], so a reader that reserves from such a count aborts the tests instead of
//! returning an error.

use std::{
    alloc::{GlobalAlloc, Layout, System},
    io, ptr,
};

use pasta::bam;

/// The largest single allocation this binary allows.
const MAX_ALLOCATION: usize = 256 << 20;

struct CappedAllocator;

// SAFETY: Every call is forwarded to `System`, or fails by returning null.
unsafe impl GlobalAlloc for CappedAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if layout.size() > MAX_ALLOCATION {
            ptr::null_mut()
        } else {
            // SAFETY: The caller upholds `GlobalAlloc::alloc`'s contract.
            unsafe { System.alloc(layout) }
        }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        if layout.size() > MAX_ALLOCATION {
            ptr::null_mut()
        } else {
            // SAFETY: The caller upholds `GlobalAlloc::alloc_zeroed`'s contract.
            unsafe { System.alloc_zeroed(layout) }
        }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: The caller upholds `GlobalAlloc::dealloc`'s contract.
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if new_size > MAX_ALLOCATION {
            ptr::null_mut()
        } else {
            // SAFETY: The caller upholds `GlobalAlloc::realloc`'s contract.
            unsafe { System.realloc(ptr, layout, new_size) }
        }
    }
}

#[global_allocator]
static ALLOCATOR: CappedAllocator = CappedAllocator;

fn put_u32_le(dst: &mut Vec<u8>, n: u32) {
    dst.extend(n.to_le_bytes());
}

/// A BAM header (uncompressed) with no text and one reference sequence, `sq0`.
fn header() -> Vec<u8> {
    let mut src = Vec::new();
    src.extend(b"BAM\x01");
    put_u32_le(&mut src, 0); // l_text
    put_u32_le(&mut src, 1); // n_ref
    put_u32_le(&mut src, 4); // ref[0].l_name
    src.extend(b"sq0\x00"); // ref[0].name
    put_u32_le(&mut src, 8); // ref[0].l_ref
    src
}

fn read_header(src: &[u8]) -> io::Result<pasta::sam::Header> {
    bam::io::Reader::from(src).read_header()
}

fn assert_unexpected_eof<T>(result: io::Result<T>) {
    match result {
        Ok(_) => panic!("expected an error"),
        Err(e) => assert_eq!(e.kind(), io::ErrorKind::UnexpectedEof, "{e}"),
    }
}

#[test]
fn test_read_header_with_huge_l_text() {
    let mut src = Vec::new();
    src.extend(b"BAM\x01");
    put_u32_le(&mut src, u32::MAX); // l_text
    src.extend(b"@HD\tVN:1.6\n");

    // The text ends at the end of the stream; then `n_ref` is missing.
    assert_unexpected_eof(read_header(&src));
}

#[test]
fn test_read_header_with_huge_n_ref() {
    let mut src = Vec::new();
    src.extend(b"BAM\x01");
    put_u32_le(&mut src, 0); // l_text
    put_u32_le(&mut src, u32::MAX); // n_ref
    put_u32_le(&mut src, 4); // ref[0].l_name
    src.extend(b"sq0\x00"); // ref[0].name
    put_u32_le(&mut src, 8); // ref[0].l_ref

    assert_unexpected_eof(read_header(&src));
}

#[test]
fn test_read_header_with_huge_l_name() {
    let mut src = Vec::new();
    src.extend(b"BAM\x01");
    put_u32_le(&mut src, 0); // l_text
    put_u32_le(&mut src, 1); // n_ref
    put_u32_le(&mut src, u32::MAX); // ref[0].l_name
    src.extend(b"sq0\x00"); // ref[0].name
    put_u32_le(&mut src, 8); // ref[0].l_ref

    assert_unexpected_eof(read_header(&src));
}

#[test]
fn test_read_record_with_huge_block_size() -> io::Result<()> {
    let mut src = header();
    put_u32_le(&mut src, u32::MAX); // block_size
    src.extend([0; 64]);

    let mut reader = bam::io::Reader::from(&src[..]);
    reader.read_header()?;

    let mut record = bam::Record::default();
    assert_unexpected_eof(reader.read_record(&mut record));

    let mut reader = bam::io::Reader::from(&src[..]);
    let header = reader.read_header()?;
    assert_unexpected_eof(reader.record_bufs(&header).next().expect("a record"));

    Ok(())
}
