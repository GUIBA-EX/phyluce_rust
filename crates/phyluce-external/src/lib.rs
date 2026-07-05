//! phyluce-external: run external bioinformatics binaries (lastz, mafft,
//! spades, ...) with a consistent, diagnosable failure mode.
//!
//! The rewrite plan calls out that legacy commands are inconsistent about
//! treating non-empty stderr as failure. This wrapper always captures both
//! streams and the exit status, and lets the caller decide the stderr policy
//! explicitly via [`StderrPolicy`] instead of silently guessing.

use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StderrPolicy {
    /// Non-empty stderr is fine (many bioinformatics tools log progress there).
    Ignore,
    /// Non-empty stderr is treated as a failure even if the exit code is 0.
    TreatAsError,
}

#[derive(Debug, thiserror::Error)]
pub enum ExternalError {
    #[error("failed to spawn '{command}': {source}")]
    Spawn {
        command: String,
        #[source]
        source: std::io::Error,
    },
    #[error("'{command}' exited with status {status}\n--- stderr ---\n{stderr}")]
    NonZeroExit {
        command: String,
        status: i32,
        stderr: String,
    },
    #[error("'{command}' wrote to stderr (treated as failure)\n--- stderr ---\n{stderr}")]
    StderrNonEmpty { command: String, stderr: String },
}

/// Full record of one external command's execution, suitable for logging.
#[derive(Debug, Clone)]
pub struct RunReport {
    pub command: String,
    pub args: Vec<String>,
    pub exit_code: Option<i32>,
    pub stdout: String,
    pub stderr: String,
    pub duration: Duration,
}

impl RunReport {
    pub fn command_line(&self) -> String {
        let mut parts = vec![self.command.clone()];
        parts.extend(self.args.iter().cloned());
        parts.join(" ")
    }
}

pub struct ExternalCommand {
    binary: PathBuf,
    args: Vec<String>,
    cwd: Option<PathBuf>,
    stderr_policy: StderrPolicy,
}

impl ExternalCommand {
    pub fn new(binary: impl Into<PathBuf>) -> Self {
        Self {
            binary: binary.into(),
            args: Vec::new(),
            cwd: None,
            stderr_policy: StderrPolicy::Ignore,
        }
    }

    pub fn arg(mut self, arg: impl Into<String>) -> Self {
        self.args.push(arg.into());
        self
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.args.extend(args.into_iter().map(Into::into));
        self
    }

    pub fn current_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.cwd = Some(dir.into());
        self
    }

    pub fn stderr_policy(mut self, policy: StderrPolicy) -> Self {
        self.stderr_policy = policy;
        self
    }

    /// Run the command to completion, capturing stdout/stderr.
    pub fn run(&self) -> Result<RunReport, ExternalError> {
        let command_str = self.binary.to_string_lossy().to_string();
        let mut cmd = Command::new(&self.binary);
        cmd.args(&self.args)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        if let Some(cwd) = &self.cwd {
            cmd.current_dir(cwd);
        }

        let start = Instant::now();
        let output = cmd.output().map_err(|source| ExternalError::Spawn {
            command: command_str.clone(),
            source,
        })?;
        let duration = start.elapsed();

        let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        let exit_code = output.status.code();

        let report = RunReport {
            command: command_str.clone(),
            args: self.args.clone(),
            exit_code,
            stdout,
            stderr: stderr.clone(),
            duration,
        };

        if !output.status.success() {
            return Err(ExternalError::NonZeroExit {
                command: report.command_line(),
                status: exit_code.unwrap_or(-1),
                stderr,
            });
        }
        if self.stderr_policy == StderrPolicy::TreatAsError && !report.stderr.trim().is_empty() {
            return Err(ExternalError::StderrNonEmpty {
                command: report.command_line(),
                stderr: report.stderr.clone(),
            });
        }

        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn captures_stdout_and_exit_code() {
        let report = ExternalCommand::new("echo").arg("hello").run().unwrap();
        assert_eq!(report.stdout.trim(), "hello");
        assert_eq!(report.exit_code, Some(0));
    }

    #[test]
    fn nonzero_exit_is_an_error() {
        let err = ExternalCommand::new("sh")
            .args(["-c", "exit 3"])
            .run()
            .unwrap_err();
        match err {
            ExternalError::NonZeroExit { status, .. } => assert_eq!(status, 3),
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn stderr_policy_can_treat_stderr_as_failure() {
        let err = ExternalCommand::new("sh")
            .args(["-c", "echo boom 1>&2"])
            .stderr_policy(StderrPolicy::TreatAsError)
            .run()
            .unwrap_err();
        assert!(matches!(err, ExternalError::StderrNonEmpty { .. }));
    }
}
