use rust_htslib::bam::{self, Read, Record, Header, HeaderView, header::HeaderRecord};
use rust_htslib::htslib;
use std::collections::BinaryHeap;
use std::cmp::Ordering;
use std::path::Path;
use anyhow;
use crate::commons::*;

// Default fields for coverage calculation (matching C++ tmerge.cpp)
pub const DEFAULT_FIELDS: u32 = htslib::sam_fields_SAM_QNAME as u32 | 
                                htslib::sam_fields_SAM_FLAG as u32 | 
                                htslib::sam_fields_SAM_RNAME as u32 | 
                                htslib::sam_fields_SAM_POS as u32 | 
                                htslib::sam_fields_SAM_CIGAR as u32 | 
                                htslib::sam_fields_SAM_AUX as u32;

#[derive(Debug, Clone)]
pub struct TBSAMReaderRecord {
    record: Record,
    file_idx: usize,
}

impl TBSAMReaderRecord {
    pub fn new(record: Record, file_idx: usize) -> Self {
        TBSAMReaderRecord { record, file_idx }
    }
    pub fn record(&self) -> &Record {
        &self.record
    }
    pub fn file_idx(&self) -> usize {
        self.file_idx
    }
}

impl PartialEq for TBSAMReaderRecord {
    fn eq(&self, other: &Self) -> bool { 
        (self.record.tid(), self.record.pos()) == (other.record.tid(), other.record.pos()) 
    }
}
impl Eq for TBSAMReaderRecord {}
impl PartialOrd for TBSAMReaderRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for TBSAMReaderRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        (other.record.tid(), other.record.pos()).cmp(&(self.record.tid(), self.record.pos()))
    }
}

pub struct TBSAMReader {
    readers: Vec<bam::Reader>,
    heap: BinaryHeap<TBSAMReaderRecord>,
    mheader: Header,
    required_fields: u32,
}

fn is_coordinate_sorted(header: &HeaderView) -> anyhow::Result<bool> {
    let header = bam::Header::from_template(header);
    if let Some(hd_vec) = header.to_hashmap().get("HD") {
        if let Some(hd_fields) = hd_vec.first() {
            if let Some(so) = hd_fields.get("SO") {
                Ok(so == "coordinate")
            } else {
                Ok(false)
            }
        }
        else {
            anyhow::bail!("SO tag not found in @HD line");
        }
    }
    else {
        anyhow::bail!("HD tag not found in @HD line");
    }
}

fn header_sq_order_eq(header1: &HeaderView, header2: &HeaderView) -> bool {
    if header1.target_count() != header2.target_count() {
        return false;
    }
    for i in 0..header1.target_count() {
        if header1.tid2name(i) != header2.tid2name(i) {
            return false;
        }
        if header1.target_len(i).unwrap() != header2.target_len(i).unwrap() {
            return false;
        }
    }
    true
}

impl TBSAMReader {
    pub fn new<P: AsRef<Path>>(files: &[P]) -> anyhow::Result<Self> {
        Self::new_with_fields(files, DEFAULT_FIELDS)
    }

    pub fn new_with_fields<P: AsRef<Path>>(files: &[P], required_fields: u32) -> anyhow::Result<Self> {
        let mut readers = Vec::new();
        let mut heap = BinaryHeap::new();
        let mut ref_header: Option<HeaderView> = None;
        let mut mheader = Header::new();

        for (i, file) in files.iter().enumerate() {
            let mut reader = bam::Reader::from_path(file)?;

            unsafe {
                htslib::hts_set_opt(
                    reader.htsfile(),
                    required_fields,
                );
            }
            
            // check header correctness
            if let Some(ref ref_hdr) = ref_header {
                let new_hdr = reader.header().clone();
                if !is_coordinate_sorted(&new_hdr)? {
                    anyhow::bail!("File {:?} is not coordinate sorted", file.as_ref());
                }
                if !header_sq_order_eq(ref_hdr, &new_hdr) {
                    anyhow::bail!("File {:?} has different reference sequence order", file.as_ref());
                }
                // merge headers
                // for every header line type (PG, RG, etc), except SQ and HD lines
                // add them to the merged header
                let new_header = Header::from_template(&new_hdr);
                merge_headers( &new_header, &mut mheader);
            } else {
                ref_header = Some(reader.header().clone());
                mheader = Header::from_template(&ref_header.clone().unwrap());
            }

            let mut record = Record::new();
            match reader.read(&mut record) {
                Some(Ok(_)) => {
                    heap.push(TBSAMReaderRecord { record, file_idx: i });
                },
                Some(Err(e)) => anyhow::bail!("Error reading file {:?}: {}", file.as_ref(), e),
                None => {},
            }
            readers.push(reader);
            // add CO line with sample data to the header
            let comment_line = format!("SAMPLE:{}", file.as_ref().canonicalize().unwrap().to_string_lossy());
            mheader.push_comment(comment_line.as_bytes());
        }

        // add PG for TieBrush
        let mut tb_pg_record = HeaderRecord::new("PG".as_bytes());
        tb_pg_record.push_tag("ID".as_bytes(), "TieBrush");
        tb_pg_record.push_tag("VN".as_bytes(), env!("CARGO_PKG_VERSION"));
        let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
        tb_pg_record.push_tag("CL".as_bytes(), cmdline);    
        mheader.push_record(&tb_pg_record);

        Ok(Self { readers, heap, mheader, required_fields })
    }

    pub fn get_header(&self) -> &Header {
        &self.mheader
    }
}
    
impl Iterator for TBSAMReader {
    type Item = TBSAMReaderRecord;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(min_rec) = self.heap.pop() {
            let file_idx = min_rec.file_idx;
            let result = min_rec.clone();

            // Try to read next record from the same file
            let reader = &mut self.readers[file_idx];
            let mut next_record = Record::new();
            match reader.read(&mut next_record) {
                Some(Ok(_)) => {
                    self.heap.push(TBSAMReaderRecord {
                        record: next_record,
                        file_idx,
                    });
                },
                Some(Err(e)) => {
                    eprintln!("Error reading file: {}", e);
                    return None;
                }
                None => {},
            }
            Some(result)
        } else {
            None
        }
    }
} 