use rust_htslib::bam::{self, Read, Record, Header, HeaderView};
use std::collections::BinaryHeap;
use std::cmp::Ordering;
use std::path::Path;
use anyhow;
use std::collections::HashMap;
use crate::commons::*;

struct SAMReaderRecord {
    record: Record,
    file_idx: usize,
}

impl PartialEq for SAMReaderRecord {
    fn eq(&self, other: &Self) -> bool { (self.record.tid(),self.record.pos()) == (other.record.tid(),other.record.pos()) }
}
impl Eq for SAMReaderRecord {}
impl PartialOrd for SAMReaderRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some((other.record.tid(),other.record.pos()).cmp(&(self.record.tid(),self.record.pos())))
    }
}
impl Ord for SAMReaderRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        (other.record.tid(),other.record.pos()).cmp(&(self.record.tid(),self.record.pos()))
    }
}

pub struct SAMReader {
    readers: Vec<bam::Reader>,
    heap: BinaryHeap<SAMReaderRecord>,
    mheader: Header,
}

fn is_coordinate_sorted(header: &HeaderView) -> anyhow::Result<bool> {
    let header = bam::Header::from_template(header);
    if let Some(hd_vec) = header.to_hashmap().get("HD") {
        if let Some(hd_fields) = hd_vec.first() {
            if let Some(sort_order) = hd_fields.get("SO") {
                Ok(true)
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

impl SAMReader {
    pub fn new<P: AsRef<Path>>(files: &[P]) -> anyhow::Result<Self> {
        let mut readers = Vec::new();
        let mut heap = BinaryHeap::new();
        let mut ref_header: Option<HeaderView> = None;

        let mut mheader = Header::new();

        for (i, file) in files.iter().enumerate() {
            let mut reader = bam::Reader::from_path(file)?;
            
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
                    heap.push(SAMReaderRecord { record, file_idx: i });
                },
                Some(Err(e)) => anyhow::bail!("Error reading file {:?}: {}", file.as_ref(), e),
                None => {},
            }
            readers.push(reader);
        }

        Ok(Self { readers, heap, mheader })
    }

    pub fn get_header(&self) -> &Header {
        &self.mheader
    }
}
    
impl Iterator for SAMReader {
    type Item = Record;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(min_rec) = self.heap.pop() {
            let file_idx = min_rec.file_idx;
            let result = min_rec.record.clone();

            // Try to read next record from the same file
            let reader = &mut self.readers[file_idx];
            let mut next_record = Record::new();
            match reader.read(&mut next_record) {
                Some(Ok(_)) => {
                    self.heap.push(SAMReaderRecord {
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