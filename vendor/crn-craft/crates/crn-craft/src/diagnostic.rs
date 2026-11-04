use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    Error,
    Warning,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Diagnostic {
    pub severity: Severity,
    pub code: String,
    /// RFC 6901 JSON Pointer. Empty string denotes the document root.
    pub path: String,
    pub message: String,
    pub hint: String,
}

impl Diagnostic {
    pub fn error(code: &str, path: &str, message: impl Into<String>, hint: &str) -> Self {
        Self {
            severity: Severity::Error,
            code: code.into(),
            path: path.into(),
            message: message.into(),
            hint: hint.into(),
        }
    }

    pub fn warning(code: &str, path: &str, message: impl Into<String>, hint: &str) -> Self {
        let mut diagnostic = Self::error(code, path, message, hint);
        diagnostic.severity = Severity::Warning;
        diagnostic
    }
}

pub fn pointer(parent: &str, key: impl std::fmt::Display) -> String {
    format!(
        "{}/{}",
        parent,
        key.to_string().replace('~', "~0").replace('/', "~1")
    )
}
