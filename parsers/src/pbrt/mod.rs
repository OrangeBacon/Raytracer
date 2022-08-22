#![allow(unused)]

mod options;
mod transform_set;

use std::{path::Path, str::Chars};

use anyhow::Result;
use geometry::{Number, Float};

use crate::pbrt::transform_set::TransformSet;

pub fn parse_pbrt(file_name: impl AsRef<Path>) -> Result<Vec<()>> {
    let mut files = vec![std::fs::read_to_string(file_name)?];

    let mut parser = Parser::<Float>::new(&mut files);

    parser.parse()
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

    fn next(&mut self) -> Option<&str> {
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
                            return Some(&tok[0..=len]);
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
                            return Some(&self.escaped_strings);
                        }
                    }
                    panic!("Non-terminated string");
                }
                '[' | ']' => {
                    return Some(&tok[0..1]);
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
                    return Some(&tok[0..len]);
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

    fn parse(&mut self) -> Result<Vec<()>> {
        Ok(vec![])
    }
}
