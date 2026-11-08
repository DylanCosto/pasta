//! MD5 (RFC 1321), for the checksums of reference sequences.
//!
//! The rounds are written so that as little as possible waits on the word computed in the step
//! before, which is what bounds MD5's speed: in the second round, the half of the mixing function
//! that does not use it is added first. This hashes about 10% faster than the `md-5` crate.

/// The digest size, in bytes.
pub(crate) const OUTPUT_SIZE: usize = 16;

const BLOCK_SIZE: usize = 64;

const INITIAL_STATE: [u32; 4] = [0x6745_2301, 0xefcd_ab89, 0x98ba_dcfe, 0x1032_5476];

/// The constants added at each step: the integer parts of `abs(sin(i + 1)) * 2^32`.
const K: [u32; 64] = [
    0xd76a_a478,
    0xe8c7_b756,
    0x2420_70db,
    0xc1bd_ceee,
    0xf57c_0faf,
    0x4787_c62a,
    0xa830_4613,
    0xfd46_9501,
    0x6980_98d8,
    0x8b44_f7af,
    0xffff_5bb1,
    0x895c_d7be,
    0x6b90_1122,
    0xfd98_7193,
    0xa679_438e,
    0x49b4_0821,
    0xf61e_2562,
    0xc040_b340,
    0x265e_5a51,
    0xe9b6_c7aa,
    0xd62f_105d,
    0x0244_1453,
    0xd8a1_e681,
    0xe7d3_fbc8,
    0x21e1_cde6,
    0xc337_07d6,
    0xf4d5_0d87,
    0x455a_14ed,
    0xa9e3_e905,
    0xfcef_a3f8,
    0x676f_02d9,
    0x8d2a_4c8a,
    0xfffa_3942,
    0x8771_f681,
    0x6d9d_6122,
    0xfde5_380c,
    0xa4be_ea44,
    0x4bde_cfa9,
    0xf6bb_4b60,
    0xbebf_bc70,
    0x289b_7ec6,
    0xeaa1_27fa,
    0xd4ef_3085,
    0x0488_1d05,
    0xd9d4_d039,
    0xe6db_99e5,
    0x1fa2_7cf8,
    0xc4ac_5665,
    0xf429_2244,
    0x432a_ff97,
    0xab94_23a7,
    0xfc93_a039,
    0x655b_59c3,
    0x8f0c_cc92,
    0xffef_f47d,
    0x8584_5dd1,
    0x6fa8_7e4f,
    0xfe2c_e6e0,
    0xa301_4314,
    0x4e08_11a1,
    0xf753_7e82,
    0xbd3a_f235,
    0x2ad7_d2bb,
    0xeb86_d391,
];

/// An MD5 hash being computed.
pub(crate) struct Md5 {
    state: [u32; 4],
    /// The bytes hashed so far.
    len: u64,
    /// The bytes of the last, incomplete block.
    buf: [u8; BLOCK_SIZE],
    buf_len: usize,
}

impl Md5 {
    pub(crate) fn new() -> Self {
        Self {
            state: INITIAL_STATE,
            len: 0,
            buf: [0; BLOCK_SIZE],
            buf_len: 0,
        }
    }

    pub(crate) fn update(&mut self, mut data: &[u8]) {
        self.len = self.len.wrapping_add(data.len() as u64);

        if self.buf_len > 0 {
            let n = (BLOCK_SIZE - self.buf_len).min(data.len());
            self.buf[self.buf_len..self.buf_len + n].copy_from_slice(&data[..n]);
            self.buf_len += n;
            data = &data[n..];

            if self.buf_len < BLOCK_SIZE {
                return;
            }

            compress(&mut self.state, &[self.buf]);
            self.buf_len = 0;
        }

        let (blocks, rest) = data.as_chunks::<BLOCK_SIZE>();
        compress(&mut self.state, blocks);

        self.buf[..rest.len()].copy_from_slice(rest);
        self.buf_len = rest.len();
    }

    pub(crate) fn finalize(mut self) -> [u8; OUTPUT_SIZE] {
        let bit_len = self.len.wrapping_mul(8);

        // A 1 bit, zeros, then the length in bits, ending a block.
        let mut block = [0; BLOCK_SIZE];
        block[..self.buf_len].copy_from_slice(&self.buf[..self.buf_len]);
        block[self.buf_len] = 0x80;

        if self.buf_len >= BLOCK_SIZE - 8 {
            compress(&mut self.state, &[block]);
            block = [0; BLOCK_SIZE];
        }

        block[BLOCK_SIZE - 8..].copy_from_slice(&bit_len.to_le_bytes());
        compress(&mut self.state, &[block]);

        let mut digest = [0; OUTPUT_SIZE];

        for (dst, word) in digest.as_chunks_mut::<4>().0.iter_mut().zip(self.state) {
            *dst = word.to_le_bytes();
        }

        digest
    }

    /// The digest of `data`.
    #[cfg(test)]
    pub(crate) fn digest(data: &[u8]) -> [u8; OUTPUT_SIZE] {
        let mut hasher = Self::new();
        hasher.update(data);
        hasher.finalize()
    }
}

fn compress(state: &mut [u32; 4], blocks: &[[u8; BLOCK_SIZE]]) {
    for block in blocks {
        let mut m = [0; 16];

        for (word, bytes) in m.iter_mut().zip(block.as_chunks::<4>().0) {
            *word = u32::from_le_bytes(*bytes);
        }

        let [mut a, mut b, mut c, mut d] = *state;

        // Each step mixes `b`, `c` and `d` into `a` with message word `m[$i]` and constant `K[$k]`,
        // rotates it by `$s`, adds `b`, and the words move along.
        macro_rules! step {
            ($early:expr, $late:expr, $i:literal, $k:literal, $s:literal) => {{
                let t = a
                    .wrapping_add(m[$i].wrapping_add(K[$k]))
                    .wrapping_add($early);
                let t = t.wrapping_add($late).rotate_left($s).wrapping_add(b);
                (a, b, c, d) = (d, t, b, c);
            }};
        }

        // F(b, c, d) = (b & c) | (!b & d)
        macro_rules! f {
            ($i:literal, $k:literal, $s:literal) => {
                step!(0, d ^ (b & (c ^ d)), $i, $k, $s)
            };
        }

        // G(b, c, d) = (b & d) | (c & !d), whose halves never overlap, so they can be added, and
        // the half without `b` first.
        macro_rules! g {
            ($i:literal, $k:literal, $s:literal) => {
                step!(c & !d, b & d, $i, $k, $s)
            };
        }

        // H(b, c, d) = b ^ c ^ d
        macro_rules! h {
            ($i:literal, $k:literal, $s:literal) => {
                step!(0, b ^ (c ^ d), $i, $k, $s)
            };
        }

        // I(b, c, d) = c ^ (b | !d)
        macro_rules! i {
            ($i:literal, $k:literal, $s:literal) => {
                step!(0, c ^ (b | !d), $i, $k, $s)
            };
        }

        f!(0, 0, 7);
        f!(1, 1, 12);
        f!(2, 2, 17);
        f!(3, 3, 22);
        f!(4, 4, 7);
        f!(5, 5, 12);
        f!(6, 6, 17);
        f!(7, 7, 22);
        f!(8, 8, 7);
        f!(9, 9, 12);
        f!(10, 10, 17);
        f!(11, 11, 22);
        f!(12, 12, 7);
        f!(13, 13, 12);
        f!(14, 14, 17);
        f!(15, 15, 22);

        g!(1, 16, 5);
        g!(6, 17, 9);
        g!(11, 18, 14);
        g!(0, 19, 20);
        g!(5, 20, 5);
        g!(10, 21, 9);
        g!(15, 22, 14);
        g!(4, 23, 20);
        g!(9, 24, 5);
        g!(14, 25, 9);
        g!(3, 26, 14);
        g!(8, 27, 20);
        g!(13, 28, 5);
        g!(2, 29, 9);
        g!(7, 30, 14);
        g!(12, 31, 20);

        h!(5, 32, 4);
        h!(8, 33, 11);
        h!(11, 34, 16);
        h!(14, 35, 23);
        h!(1, 36, 4);
        h!(4, 37, 11);
        h!(7, 38, 16);
        h!(10, 39, 23);
        h!(13, 40, 4);
        h!(0, 41, 11);
        h!(3, 42, 16);
        h!(6, 43, 23);
        h!(9, 44, 4);
        h!(12, 45, 11);
        h!(15, 46, 16);
        h!(2, 47, 23);

        i!(0, 48, 6);
        i!(7, 49, 10);
        i!(14, 50, 15);
        i!(5, 51, 21);
        i!(12, 52, 6);
        i!(3, 53, 10);
        i!(10, 54, 15);
        i!(1, 55, 21);
        i!(8, 56, 6);
        i!(15, 57, 10);
        i!(6, 58, 15);
        i!(13, 59, 21);
        i!(4, 60, 6);
        i!(11, 61, 10);
        i!(2, 62, 15);
        i!(9, 63, 21);

        for (s, x) in state.iter_mut().zip([a, b, c, d]) {
            *s = s.wrapping_add(x);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(digest: [u8; OUTPUT_SIZE]) -> String {
        digest.iter().map(|b| format!("{b:02x}")).collect()
    }

    #[test]
    fn test_digest() {
        // RFC 1321, A.5 "Test suite".
        for (src, expected) in [
            (&b""[..], "d41d8cd98f00b204e9800998ecf8427e"),
            (b"a", "0cc175b9c0f1b6a831c399e269772661"),
            (b"abc", "900150983cd24fb0d6963f7d28e17f72"),
            (b"message digest", "f96b697d7cb7938d525a2f31aaf161d0"),
            (
                b"abcdefghijklmnopqrstuvwxyz",
                "c3fcd3d76192e4007dfb496cca67e13b",
            ),
            (
                b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
                "d174ab98d277d9f5a5611c2c9f419d9f",
            ),
            (
                b"12345678901234567890123456789012345678901234567890123456789012345678901234567890",
                "57edf4a22be3c955ac49da2e2107b67a",
            ),
        ] {
            assert_eq!(hex(Md5::digest(src)), expected, "{src:?}");
        }

        // A million `a`s.
        assert_eq!(
            hex(Md5::digest(&[b'a'; 1_000_000])),
            "7707d6ae4e027c70eea2a935c2296f21"
        );
    }

    #[test]
    fn test_update_in_parts() {
        let src: Vec<u8> = (0..1000u32).map(|i| (i * 7 + i / 13) as u8).collect();
        let expected = Md5::digest(&src);

        // Parts ending in and across blocks, and lengths near the padding boundary.
        for part_len in [1, 3, 55, 56, 63, 64, 65, 127, 999] {
            let mut hasher = Md5::new();

            for part in src.chunks(part_len) {
                hasher.update(part);
            }

            assert_eq!(hasher.finalize(), expected, "part length {part_len}");
        }

        for len in 50..=70 {
            let mut hasher = Md5::new();
            hasher.update(&src[..len]);
            hasher.update(&[]);
            assert_eq!(hasher.finalize(), Md5::digest(&src[..len]), "length {len}");
        }
    }
}
