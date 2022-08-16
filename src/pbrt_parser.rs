use std::{path::Path, str::Chars};

use anyhow::Result;

pub fn parse_pbrt(file_name: impl AsRef<Path>) -> Result<Vec<Commands>> {
    let mut files = vec![std::fs::read_to_string(file_name)?];

    let mut parser = Parser {
        files: &mut files,
    };

    Ok(parser.parse()?)
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
                    while let Some(ch) = self.chars.next() {
                        if ch == '\n' || ch == '\r' {
                            break;
                        }
                    }
                }
                '"' => {
                    let mut len = 0;
                    while let Some(ch) = self.chars.next() {
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

#[derive(Debug)]
pub enum Commands {
    AttributeBegin,
    AttributeEnd,
    ActiveTransform,
    AreaLightSource,
    Accelerator,
    ConcatTransform,
    CoordinateSystem,
    CoordSysTransform,
    Camera,
    Film,
    Integrator,
    Include,
    Identity,
    LightSource,
    LookAt,
    MakeNamedMaterial,
    MakeNamedMedium,
    Material,
    MediumInterface,
    NamedMaterial,
    ObjectBegin,
    ObjectEnd,
    ObjectInstance,
    PixelFilter,
    ReverseOrientation,
    Shape,
    Sampler,
    Scale,
    TransformBegin,
    TransformEnd,
    Transform,
    Translate,
    TransformTimes,
    Texture,
    WorldBegin,
    WorldEnd,
}

struct Parser<'a> {
    files: &'a mut [String],
}

impl Parser<'_> {
    fn parse(&mut self) -> Result<Vec<Commands>> {
        Ok(vec![])
    }
}