#![allow(unused)]

mod ast;
mod options;
mod transform_set;

pub use ast::*;

use std::{fmt::Display, path::Path, str::Chars};

use geometry::{Float, Number};

use self::transform_set::TransformSet;

/// Error messages encountered while parsing a pbrt file
#[derive(Debug)]
pub enum PbrtError {
    IOError(std::io::Error),
    NonTerminatedString,
}

impl From<std::io::Error> for PbrtError {
    fn from(err: std::io::Error) -> Self {
        Self::IOError(err)
    }
}

impl Display for PbrtError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PbrtError::IOError(err) => {
                f.write_fmt(format_args!("IO Error while parsing pbrt file: {}", err))
            }
            PbrtError::NonTerminatedString => f.write_str("Non-terminated string"),
        }
    }
}

impl std::error::Error for PbrtError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            PbrtError::IOError(err) => Some(err),
            _ => None,
        }
    }
}

impl<T: Number> PbrtFile<T> {
    /// Parse a pbrt file
    pub fn parse(file_name: impl AsRef<Path>) -> Result<PbrtFile<T>, PbrtError> {
        let mut files = vec![std::fs::read_to_string(file_name)?];

        let mut parser = Parser::new(&mut files);

        parser.parse()
    }
}

struct Tokeniser<'a> {
    chars: Chars<'a>,
    escaped_strings: String,
}

impl<'a> Tokeniser<'a> {
    fn new(data: &'a str) -> Self {
        Self {
            chars: data.chars(),
            escaped_strings: String::new(),
        }
    }

    fn next(&mut self) -> Option<Result<&str, PbrtError>> {
        loop {
            let tok = self.chars.as_str();
            let next = self.chars.next()?;

            match next {
                ' ' | '\n' | '\t' | '\r' => {}
                '#' => {
                    for ch in self.chars.by_ref() {
                        if ch == '\n' || ch == '\r' {
                            break;
                        }
                    }
                }
                '"' => {
                    let mut len = 0;
                    for ch in self.chars.by_ref() {
                        len += ch.len_utf8();

                        if ch == '"' {
                            return Some(Ok(&tok[0..=len]));
                        }
                        if ch == '\\' {
                            break;
                        }
                    }
                    if tok.len() <= len {
                        panic!("Non-terminated string");
                    }
                    self.escaped_strings.clear();
                    self.escaped_strings.push_str(&tok[0..len]);
                    while let Some(ch) = self.chars.next() {
                        if ch == '\\' {
                            let ch = match self.chars.next() {
                                Some('b') => '\x08',
                                Some('f') => '\x0C',
                                Some('n') => '\n',
                                Some('r') => '\r',
                                Some('t') => '\t',
                                Some('\\') => '\\',
                                Some('\'') => '\'',
                                Some('\"') => '\"',
                                Some(_) => panic!("Unexpected escape sequence"),
                                None => panic!("Unexpected end of file inside string escape"),
                            };

                            self.escaped_strings.push(ch);
                        } else {
                            self.escaped_strings.push(ch);
                        }

                        if ch == '"' {
                            return Some(Ok(&self.escaped_strings));
                        }
                    }
                    panic!("Non-terminated string");
                }
                '[' | ']' => {
                    return Some(Ok(&tok[0..1]));
                }
                _ => {
                    let mut len = next.len_utf8();
                    for ch in self.chars.clone() {
                        if " \n\t\r\"[]".contains(ch) {
                            break;
                        }
                        len += ch.len_utf8();
                        self.chars.next();
                    }
                    return Some(Ok(&tok[0..len]));
                }
            }
        }
    }
}

struct Parser<'a, T: Number> {
    files: &'a mut [String],
    transforms: TransformSet<T>,
}

impl<'a, T: Number> Parser<'a, T> {
    fn new(files: &'a mut [String]) -> Self {
        Self {
            files,
            transforms: TransformSet::new(),
        }
    }

    fn parse(&mut self) -> Result<PbrtFile<T>, PbrtError> {
        Ok(PbrtFile::default())
    }
}
