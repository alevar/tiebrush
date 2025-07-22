use rust_htslib::bam::{self, Read, Record};
use std::collections::BinaryHeap;
use std::cmp::Ordering;
use std::path::Path;
use anyhow;

struct SAMReaderRecord {
    record: Record,
    file_idx: usize,
}

impl PartialEq for SAMReaderRecord {
    fn eq(&self, other: &Self) -> bool { self.record.pos() == other.record.pos() }
}
impl Eq for SAMReaderRecord {}
impl PartialOrd for SAMReaderRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(other.record.pos().cmp(&self.record.pos()))
    }
}
impl Ord for SAMReaderRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        other.record.pos().cmp(&self.record.pos())
    }
}

pub struct SAMReader {
    readers: Vec<bam::Reader>,
    heap: BinaryHeap<SAMReaderRecord>,
}

fn is_coordinate_sorted(reader: &bam::Reader) -> anyhow::Result<bool> {
    let header = bam::Header::from_template(reader.header());
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

impl SAMReader {
    pub fn new<P: AsRef<Path>>(files: &[P]) -> anyhow::Result<Self> {
        let mut readers = Vec::new();
        let mut heap = BinaryHeap::new();

        for (i, file) in files.iter().enumerate() {
            let mut reader = bam::Reader::from_path(file)?;
            if !is_coordinate_sorted(&reader)? {
                anyhow::bail!("File {:?} is not coordinate sorted", file.as_ref());
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

        Ok(Self { readers, heap })
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