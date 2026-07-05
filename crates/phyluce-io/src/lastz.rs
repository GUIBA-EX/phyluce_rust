//! Reader for LASTZ `--format=general-` output, mirroring `phyluce/lastz.py`'s
//! `Reader` class field-for-field (including its `>`-stripping of name1/name2
//! and int/percent coercion of specific columns).

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum LastzError {
    #[error("{0}")]
    Io(#[from] io::Error),
    #[error("line {line}: expected {expected} tab-separated fields, found {found}")]
    FieldCount {
        line: usize,
        expected: usize,
        found: usize,
    },
    #[error("line {line}: could not parse integer field {field}: {value:?}")]
    BadInt {
        line: usize,
        field: usize,
        value: String,
    },
    #[error("line {line}: could not parse percent field {field}: {value:?}")]
    BadPercent {
        line: usize,
        field: usize,
        value: String,
    },
}

/// One row of `lastz --format=general-:score,name1,strand1,zstart1,end1,
/// length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,
/// continuity[,coverage]` output.
#[derive(Debug, Clone, PartialEq)]
pub struct LastzMatch {
    pub score: String,
    pub name1: String,
    pub strand1: String,
    pub zstart1: i64,
    pub end1: i64,
    pub length1: i64,
    pub name2: String,
    pub strand2: String,
    pub zstart2: i64,
    pub end2: i64,
    pub length2: i64,
    pub diff: String,
    pub cigar: String,
    pub identity: String,
    pub percent_identity: f64,
    pub continuity: String,
    pub percent_continuity: f64,
    /// Only present when reading long-format LASTZ output.
    pub coverage: Option<String>,
    pub percent_coverage: Option<f64>,
}

pub fn read_lastz(path: &Path, long_format: bool) -> Result<Vec<LastzMatch>, LastzError> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    let expected = if long_format { 19 } else { 17 };
    let mut matches = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let lineno = i + 1;
        let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        if fields.len() != expected {
            return Err(LastzError::FieldCount {
                line: lineno,
                expected,
                found: fields.len(),
            });
        }

        let get_int = |fields: &[String], idx: usize| -> Result<i64, LastzError> {
            fields[idx].parse::<i64>().map_err(|_| LastzError::BadInt {
                line: lineno,
                field: idx,
                value: fields[idx].clone(),
            })
        };
        let get_percent = |fields: &[String], idx: usize| -> Result<f64, LastzError> {
            fields[idx]
                .trim_end_matches('%')
                .parse::<f64>()
                .map_err(|_| LastzError::BadPercent {
                    line: lineno,
                    field: idx,
                    value: fields[idx].clone(),
                })
        };

        fields[1] = fields[1].trim_start_matches('>').to_string();
        fields[6] = fields[6].trim_start_matches('>').to_string();

        let m = LastzMatch {
            score: fields[0].clone(),
            name1: fields[1].clone(),
            strand1: fields[2].clone(),
            zstart1: get_int(&fields, 3)?,
            end1: get_int(&fields, 4)?,
            length1: get_int(&fields, 5)?,
            name2: fields[6].clone(),
            strand2: fields[7].clone(),
            zstart2: get_int(&fields, 8)?,
            end2: get_int(&fields, 9)?,
            length2: get_int(&fields, 10)?,
            diff: fields[11].clone(),
            cigar: fields[12].clone(),
            identity: fields[13].clone(),
            percent_identity: get_percent(&fields, 14)?,
            continuity: fields[15].clone(),
            percent_continuity: get_percent(&fields, 16)?,
            coverage: long_format.then(|| fields[17].clone()),
            percent_coverage: if long_format {
                Some(get_percent(&fields, 18)?)
            } else {
                None
            },
        };
        matches.push(m);
    }
    Ok(matches)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp(name: &str, contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("phyluce-io-lastz-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    #[test]
    fn parses_short_format_row() {
        let line = "11415\t>NODE_1_length_1151_cov_27.333029\t+\t576\t696\t120\t\
                     >uce-553_p1 |source:faircloth,probes-id:697,probes-locus:553,probes-probe:1\t\
                     -\t0\t120\t120\t\
                     ........................................................................................................................\t\
                     120M\t120/120\t100.0%\t120/120\t100.0%\n";
        let path = write_temp("one.lastz", line);
        let matches = read_lastz(&path, false).unwrap();
        assert_eq!(matches.len(), 1);
        let m = &matches[0];
        assert_eq!(m.name1, "NODE_1_length_1151_cov_27.333029");
        assert_eq!(m.zstart1, 576);
        assert_eq!(m.end1, 696);
        assert_eq!(
            m.name2,
            "uce-553_p1 |source:faircloth,probes-id:697,probes-locus:553,probes-probe:1"
        );
        assert_eq!(m.percent_identity, 100.0);
        assert_eq!(m.percent_continuity, 100.0);
        assert!(m.coverage.is_none());
    }
}
