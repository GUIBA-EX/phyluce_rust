//! Helpers for SQL fragments that cannot be bound as parameters.
//!
//! rusqlite parameters cover values, but table and column names still have to
//! be interpolated into SQL text. Always quote those identifiers here instead
//! of formatting raw names into statements.

pub fn ident(name: &str) -> String {
    format!("\"{}\"", name.replace('"', "\"\""))
}

pub fn qualified_ident(name: &str) -> String {
    name.split('.').map(ident).collect::<Vec<_>>().join(".")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quotes_identifiers() {
        assert_eq!(ident("taxon-name"), "\"taxon-name\"");
        assert_eq!(ident("a\"b"), "\"a\"\"b\"");
    }

    #[test]
    fn quotes_qualified_identifiers_part_by_part() {
        assert_eq!(
            qualified_ident("extended.match_map"),
            "\"extended\".\"match_map\""
        );
    }
}
