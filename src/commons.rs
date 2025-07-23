use rust_htslib::bam::{Record,record::{Aux}, Header, HeaderView, header::HeaderRecord};
use std::fmt::{self, Display};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

// convert strand to char "+" or "-" or "."
impl Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::Unknown => ".",
        })
    }
}

pub fn flags_set(record: &Record, flags: u16) -> bool {
    (record.flags() & flags) == flags
}

pub fn get_yc_tag(record: &Record) -> anyhow::Result<Option<i32>> {
    match record.aux(b"YC") {
        Ok(yc_val) => match yc_val {
            Aux::I8(val) => Ok(Some(val as i32)),
            Aux::I16(val) => Ok(Some(val as i32)),
            Aux::I32(val) => Ok(Some(val as i32)),
            Aux::U8(val) => Ok(Some(val as i32)),
            Aux::U16(val) => Ok(Some(val as i32)),
            Aux::U32(val) => Ok(Some(val as i32)),
            _ => anyhow::bail!("Value for YC tag is not an integer."),
        },
        _ => Ok(None),
    }
}

pub fn merge_headers(new_header: &Header, mheader: &mut Header) {
    for (hdr_key, hdr_vec) in new_header.to_hashmap() {
        match hdr_key.as_str() {
            "HD" | "SQ" => continue,
            _ => {
                for map in &hdr_vec {
                    let mut hdr_record = HeaderRecord::new(hdr_key.as_bytes());
                    for (key, value) in map.iter() {
                        hdr_record.push_tag(key.as_bytes(), value);
                    }
                    mheader.push_record(&hdr_record);
                }
            }
        }
    }
    for comment in new_header.comments() {
        mheader.push_comment(comment.as_bytes());
    }
}

pub fn get_strand(record: &Record) -> anyhow::Result<Strand> {
    // Try XS tag first
    match record.aux(b"XS") {
        Ok(xs_val) => match xs_val {
            Aux::Char(c) => {
                match c {
                    b'+' => {return Ok(Strand::Plus)},
                    b'-' => {return Ok(Strand::Minus)},
                    _ => {return Ok(Strand::Unknown)},
                }
            },
            _ => {anyhow::bail!("XS tag is not a character.")},
        },
        _ => {},
    }

    // Try minimap2's "ts" tag
    let is_reverse = record.is_reverse();
    match record.aux(b"ts") {
        Ok(ts_val) => match ts_val {
            Aux::Char(c) => {
                match c {
                    b'+' => {
                        if is_reverse {
                            return Ok(Strand::Minus)
                        }
                        else {
                            return Ok(Strand::Plus)
                        }
                    },
                    b'-' => {
                        if is_reverse {
                            return Ok(Strand::Plus)
                        }
                        else {
                            return Ok(Strand::Minus)
                        }
                    },
                    _ => {return Ok(Strand::Unknown)},
                }
            },
            _ => {anyhow::bail!("ts tag is not a character.")},
        },
        _ => {},
    }

    // Default: unknown
    Ok(Strand::Unknown)
}