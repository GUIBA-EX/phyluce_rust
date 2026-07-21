//! CLI wiring for `phyluce genetrees get-mean-bootrep-support`, mirroring
//! `phyluce_genetrees_get_mean_bootrep_support`.

use std::path::Path;

use anyhow::Context;
use phyluce_genetrees::newick::{parse, Node};

use crate::conf::parse_ini;

fn collect_support_values(node: &Node, out: &mut Vec<f64>) {
    if !node.is_leaf() {
        if let Some(label) = &node.label {
            if let Ok(v) = label.parse::<f64>() {
                out.push(v);
            }
        }
        for c in &node.children {
            collect_support_values(c, out);
        }
    }
}

fn support_values(tree: &Node, locus: &str) -> anyhow::Result<Vec<f64>> {
    let mut support = Vec::new();
    collect_support_values(tree, &mut support);
    anyhow::ensure!(
        !support.is_empty(),
        "tree for locus {locus:?} has no numeric bootstrap support values"
    );
    Ok(support)
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

/// Sample standard deviation (ddof=1) based 95% CI, matching
/// `1.96 * (std(ddof=1) / sqrt(n))`.
fn ci95(values: &[f64]) -> f64 {
    if values.len() <= 1 {
        return f64::NAN;
    }
    let m = mean(values);
    let variance =
        values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (values.len() as f64 - 1.0);
    1.96 * (variance.sqrt() / (values.len() as f64).sqrt())
}

/// Mirrors the legacy script's hardcoded `open("outfile.csv", "w")` --
/// there is no `--output` CLI argument in the Python original.
pub fn run(trees_dir: &Path, config: &Path) -> anyhow::Result<()> {
    let conf_text = std::fs::read_to_string(config)
        .with_context(|| format!("reading config file {}", config.display()))?;
    let sections = parse_ini(&conf_text);

    let mut out = String::from("set,mean bootrep support\n");
    for (section, entries) in &sections {
        let mut section_means = Vec::new();
        for (locus, _) in entries {
            let tree_path = trees_dir.join(locus).join("RAxML_bipartitions.FINAL");
            let text = std::fs::read_to_string(&tree_path)
                .with_context(|| format!("reading bootstrap tree file {}", tree_path.display()))?;
            let tree = parse(&text)?;
            let support = support_values(&tree, locus)?;
            section_means.push(mean(&support));
        }
        anyhow::ensure!(
            !section_means.is_empty(),
            "section {section:?} contains no loci"
        );
        let section_mean = mean(&section_means);
        let ci = ci95(&section_means);
        crate::cli_info!("{section},{},{section_mean},{ci}", section_means.len());
        for value in &section_means {
            out.push_str(&format!("{section},{value}\n"));
        }
    }
    std::fs::write("outfile.csv", out).context("writing bootstrap support summary outfile.csv")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_trees_without_numeric_support() {
        let unsupported = parse("((a,b),(c,d));").unwrap();
        assert!(support_values(&unsupported, "uce-1").is_err());

        let supported = parse("((a,b)95,(c,d)87);").unwrap();
        assert_eq!(support_values(&supported, "uce-1").unwrap(), [95.0, 87.0]);
    }
}
