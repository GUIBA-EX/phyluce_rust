//! Locus/taxon splitting mirroring `phyluce_assembly_explode_get_fastas_file`.

/// Mirrors the default (`--by-taxon` unset) grouping key: the first
/// `split_char`-delimited token of the record id (the locus name).
pub fn locus_key(id: &str, split_char: &str) -> String {
    id.split(split_char).next().unwrap_or(id).to_string()
}

/// Mirrors `--by-taxon`'s grouping key: every token after the first,
/// rejoined with `-` (organism names use `_` internally, so this
/// reconstitutes e.g. `alligator_mississippiensis` as
/// `alligator-mississippiensis`).
pub fn taxon_key(id: &str, split_char: &str) -> String {
    let parts: Vec<&str> = id.split(split_char).collect();
    if parts.len() <= 1 {
        String::new()
    } else {
        parts[1..].join("-")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locus_key_takes_first_token() {
        assert_eq!(
            locus_key("uce-1732_alligator_mississippiensis", "_"),
            "uce-1732"
        );
    }

    #[test]
    fn taxon_key_rejoins_remaining_tokens_with_dash() {
        assert_eq!(
            taxon_key("uce-1732_alligator_mississippiensis", "_"),
            "alligator-mississippiensis"
        );
    }
}
