use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use crate::fasta;

use super::{Fields, MdNm, Reader, Tags};

/// A CRAM reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
    md_nm: MdNm,
    required_fields: Fields,
    required_tags: Tags,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::reader::Builder;
    /// use pasta::fasta;
    ///
    /// let reference_sequence_repository = fasta::Repository::default();
    ///
    /// let builder = Builder::default()
    ///     .set_reference_sequence_repository(reference_sequence_repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets whether decoded records get the MD and NM tags the file left out.
    ///
    /// By default, they are regenerated, as HTSlib does.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::reader::{Builder, MdNm};
    /// let builder = Builder::default().set_md_nm(MdNm::AsStored);
    /// ```
    pub fn set_md_nm(mut self, md_nm: MdNm) -> Self {
        self.md_nm = md_nm;
        self
    }

    /// Sets the record fields the reader decodes; see [`Fields`].
    ///
    /// By default, every field is decoded.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::reader::{Builder, Fields};
    /// let builder = Builder::default().set_required_fields(Fields::FLAGS);
    /// ```
    pub fn set_required_fields(mut self, fields: Fields) -> Self {
        self.required_fields = fields;
        self
    }

    /// Sets the tags the reader decodes with [`Fields::DATA`]; see [`Tags`].
    ///
    /// By default, every tag is decoded.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::reader::{Builder, Tags};
    /// use pasta::sam::alignment::record::data::field::Tag;
    ///
    /// let builder = Builder::default().set_required_tags(Tags::Only(vec![Tag::EDIT_DISTANCE]));
    /// ```
    pub fn set_required_tags(mut self, tags: Tags) -> Self {
        self.required_tags = tags;
        self
    }

    /// Builds a CRAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use pasta::cram::io::reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.cram")?;
    /// # Ok::<_, std::io::Error>(())
    /// ````
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        File::open(src).map(|file| self.build_from_reader(file))
    }

    /// Builds a CRAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram::io::reader::Builder;
    /// let reader = Builder::default().build_from_reader(io::empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: Read,
    {
        Reader {
            inner: reader,
            reference_sequence_repository: self.reference_sequence_repository,
            md_nm: self.md_nm,
            required_fields: self.required_fields,
            required_tags: self.required_tags,
        }
    }
}
