use rust_htslib::bam::{Record,record::{Aux}, Header, HeaderView, header::HeaderRecord};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
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