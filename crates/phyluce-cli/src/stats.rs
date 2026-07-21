//! Summary statistics matching `phyluce_assembly_get_fasta_lengths`'s
//! numpy-based report (mean, sample stderr, min/max/median, contigs > 1kb).

pub struct LengthReport {
    pub count: usize,
    pub sum: u64,
    pub avg: f64,
    pub stderr: f64,
    /// True when `count <= 1`: the legacy Python sets `std_error = 0` (an
    /// `int`) in that branch, which formats as `"0"`, not `"0.0"`.
    pub stderr_is_int_zero: bool,
    pub min: Option<u64>,
    pub max: Option<u64>,
    pub median: f64,
    pub ge_1000: usize,
}

impl LengthReport {
    pub fn from_lengths(lengths: &[usize]) -> Self {
        let count = lengths.len();
        let sum: u64 = lengths.iter().map(|&x| x as u64).sum();
        let avg = if count > 0 {
            sum as f64 / count as f64
        } else {
            f64::NAN
        };
        let stderr = if count > 1 {
            let mean = avg;
            let variance: f64 = lengths
                .iter()
                .map(|&x| {
                    let d = x as f64 - mean;
                    d * d
                })
                .sum::<f64>()
                / (count as f64 - 1.0);
            variance.sqrt() / (count as f64).sqrt()
        } else {
            0.0
        };
        let min = lengths.iter().min().map(|&x| x as u64);
        let max = lengths.iter().max().map(|&x| x as u64);
        let median = median_of(lengths);
        let ge_1000 = lengths.iter().filter(|&&x| x >= 1000).count();

        Self {
            count,
            sum,
            avg,
            stderr,
            stderr_is_int_zero: count <= 1,
            min,
            max,
            median,
            ge_1000,
        }
    }

    fn fmt_stderr(&self) -> String {
        if self.stderr_is_int_zero {
            "0".to_string()
        } else {
            fmt_float(self.stderr)
        }
    }

    pub fn to_human_report(&self) -> String {
        format!(
            "Reads:\t\t{}\nBp:\t\t{}\nAvg. len:\t{}\nSTDERR len:\t{}\nMin. len:\t{}\nMax. len:\t{}\nMedian len:\t{}\nContigs > 1kb:\t{}\n",
            group_thousands(&self.count.to_string()),
            group_thousands(&self.sum.to_string()),
            group_thousands(&fmt_float(self.avg)),
            group_thousands(&self.fmt_stderr()),
            group_thousands(&self.min.map(|v| v.to_string()).unwrap_or_default()),
            group_thousands(&self.max.map(|v| v.to_string()).unwrap_or_default()),
            group_thousands(&fmt_float(self.median)),
            group_thousands(&self.ge_1000.to_string()),
        )
    }

    /// Mirrors the legacy CSV branch, including its `Div/0` fallback (a
    /// shortened row) when the input is empty and computing min/max panics
    /// in the Python original.
    pub fn to_csv_row(&self, basename: &str) -> String {
        if self.count == 0 {
            format!(
                "{},{},Div/0,Div/0,Div/0,Div/0,{}",
                basename, self.count, self.ge_1000
            )
        } else {
            format!(
                "{},{},{},{},{},{},{},{},{}",
                basename,
                self.count,
                self.sum,
                fmt_float(self.avg),
                self.fmt_stderr(),
                self.min.unwrap(),
                self.max.unwrap(),
                fmt_float(self.median),
                self.ge_1000
            )
        }
    }
}

/// Same shape of stats as [`LengthReport`], but matching
/// `phyluce_assembly_get_fastq_lengths`'s report: no "Contigs > 1kb" line
/// (that script never computes it), a literal `"All files in dir with "`
/// CSV prefix, and no count<=1 zero-stderr special case (the legacy script
/// always divides by `sqrt(len(lengths))` unconditionally; real FASTQ
/// directories never hit the single-read edge case).
pub struct FastqLengthReport {
    pub count: usize,
    pub sum: u64,
    pub avg: f64,
    pub stderr: f64,
    pub min: Option<u64>,
    pub max: Option<u64>,
    pub median: f64,
}

impl FastqLengthReport {
    pub fn from_lengths(lengths: &[usize]) -> Self {
        let count = lengths.len();
        let sum: u64 = lengths.iter().map(|&x| x as u64).sum();
        let avg = if count > 0 {
            sum as f64 / count as f64
        } else {
            f64::NAN
        };
        let variance: f64 = lengths
            .iter()
            .map(|&x| {
                let d = x as f64 - avg;
                d * d
            })
            .sum::<f64>()
            / (count as f64 - 1.0);
        let stderr = variance.sqrt() / (count as f64).sqrt();
        Self {
            count,
            sum,
            avg,
            stderr,
            min: lengths.iter().min().map(|&x| x as u64),
            max: lengths.iter().max().map(|&x| x as u64),
            median: median_of(lengths),
        }
    }

    pub fn to_human_report(&self) -> String {
        if self.count == 0 {
            // Matches `to_csv_row`'s `Div/0` handling for the same input:
            // without this, avg/stderr/median compute as literal `NaN`
            // here while the CSV path already has a clean fallback for an
            // empty/all-empty FASTQ directory.
            return "Reads:\t\t0\nBp:\t\t0\nAvg. len:\tDiv/0\nSTDERR len:\tDiv/0\nMin. len:\tDiv/0\nMax. len:\tDiv/0\nMedian len:\tDiv/0\n".to_string();
        }
        format!(
            "Reads:\t\t{}\nBp:\t\t{}\nAvg. len:\t{}\nSTDERR len:\t{}\nMin. len:\t{}\nMax. len:\t{}\nMedian len:\t{}\n",
            group_thousands(&self.count.to_string()),
            group_thousands(&self.sum.to_string()),
            group_thousands(&fmt_float(self.avg)),
            group_thousands(&fmt_float(self.stderr)),
            group_thousands(&self.min.map(|v| v.to_string()).unwrap_or_default()),
            group_thousands(&self.max.map(|v| v.to_string()).unwrap_or_default()),
            group_thousands(&fmt_float(self.median)),
        )
    }

    /// Mirrors the legacy CSV branch's literal `"All files in dir with "`
    /// prefix and its `Div/0` fallback for an empty/degenerate input.
    pub fn to_csv_row(&self, basename: &str) -> String {
        if self.count == 0 {
            format!(
                "All files in dir with {},{},Div/0,Div/0,Div/0,Div/0",
                basename, self.count
            )
        } else {
            format!(
                "All files in dir with {},{},{},{},{},{},{},{}",
                basename,
                self.count,
                self.sum,
                fmt_float(self.avg),
                fmt_float(self.stderr),
                self.min.unwrap(),
                self.max.unwrap(),
                fmt_float(self.median),
            )
        }
    }
}

fn median_of(lengths: &[usize]) -> f64 {
    if lengths.is_empty() {
        return f64::NAN;
    }
    let mut sorted: Vec<usize> = lengths.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2] as f64
    } else {
        (sorted[n / 2 - 1] as f64 + sorted[n / 2] as f64) / 2.0
    }
}

/// Format a float the way Python's `str()`/`repr()` would for use inside a
/// `"{:,}"`-style grouped string: shortest round-trip representation, no
/// forced trailing zeros for whole numbers (Python prints `123.0`, not `123`).
fn fmt_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x == x.trunc() {
        format!("{:.1}", x)
    } else {
        format!("{}", x)
    }
}

/// Insert `,` every 3 digits in the integer part of a numeric string,
/// preserving a leading `-` and any fractional part.
fn group_thousands(s: &str) -> String {
    let (sign, rest) = match s.strip_prefix('-') {
        Some(r) => ("-", r),
        None => ("", s),
    };
    let (int_part, frac_part) = match rest.split_once('.') {
        Some((i, f)) => (i, Some(f)),
        None => (rest, None),
    };

    let bytes = int_part.as_bytes();
    let mut grouped = String::new();
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            grouped.push(',');
        }
        grouped.push(*b as char);
    }

    match frac_part {
        Some(f) => format!("{sign}{grouped}.{f}"),
        None => format!("{sign}{grouped}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_thousands_correctly() {
        assert_eq!(group_thousands("1234567"), "1,234,567");
        assert_eq!(group_thousands("123"), "123");
        assert_eq!(group_thousands("1234.5678"), "1,234.5678");
        assert_eq!(group_thousands("-1234"), "-1,234");
    }

    #[test]
    fn median_odd_and_even() {
        assert_eq!(median_of(&[1, 2, 3]), 2.0);
        assert_eq!(median_of(&[1, 2, 3, 4]), 2.5);
    }

    #[test]
    fn fastq_human_report_shows_div0_not_nan_for_empty_input() {
        let report = FastqLengthReport::from_lengths(&[]);
        let human = report.to_human_report();
        assert!(
            !human.contains("nan"),
            "human report contained 'nan': {human}"
        );
        assert!(
            human.contains("Div/0"),
            "human report missing Div/0: {human}"
        );
    }

    #[test]
    fn report_basic_stats() {
        let report = LengthReport::from_lengths(&[100, 200, 300, 1200]);
        assert_eq!(report.count, 4);
        assert_eq!(report.sum, 1800);
        assert_eq!(report.min, Some(100));
        assert_eq!(report.max, Some(1200));
        assert_eq!(report.ge_1000, 1);
    }
}
