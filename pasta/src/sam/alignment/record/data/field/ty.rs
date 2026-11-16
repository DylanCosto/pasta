/// An alignment record data field value type.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Type {
    /// Character (`A`).
    Character,
    /// 8-bit integer (`c`).
    Int8,
    /// 8-bit unsigned integer (`C`).
    UInt8,
    /// 16-bit integer (`s`).
    Int16,
    /// 16-bit unsigned integer (`S`).
    UInt16,
    /// 32-bit integer (`i`).
    Int32,
    /// 32-bit unsigned integer (`I`).
    UInt32,
    /// Single-precision floating-point (`f`).
    Float,
    /// Double-precision floating-point (`d`).
    ///
    /// Not in the SAM specification, but HTSlib reads and writes it, and samtools writes it to
    /// BAM.
    Double,
    /// String (`Z`).
    String,
    /// Hex string (`H`).
    Hex,
    /// Array (`B`).
    Array,
}
