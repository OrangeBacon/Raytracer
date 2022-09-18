use std::str::Chars;

/// A single token
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Token<'a> {
    Def,
    Extern,
    Identifier(&'a str),
    Number(f64),
    Char(char),
}

/// A 'not iterator' to create tokens from a string
pub struct Scanner<'a> {
    chars: Chars<'a>,
}

impl<'a> Scanner<'a> {
    /// Create a new scanner from the provided source code
    pub fn new(source: &'a str) -> Self {
        Self {
            chars: source.chars(),
        }
    }

    /// Print all the tokens in the source string to stdout
    pub fn print_tokens(mut self) {
        while let Some(tok) = self.next() {
            println!("{:?}", tok);
        }
    }

    /// Get the next token from the source code
    pub fn next(&mut self) -> Option<Token<'a>> {
        loop {
            let tok = self.chars.as_str();
            let ch = self.chars.next()?;

            break match ch {
                '#' => {
                    // skip comments
                    while !matches!(self.chars.next(), Some('\n' | '\r')) {}
                    continue;
                }
                '0'..='9' => self.parse_number(tok),
                _ if ch.is_whitespace() => continue,
                _ if ch.is_alphabetic() => self.parse_identifier(tok, ch),
                _ => Some(Token::Char(ch)),
            };
        }
    }

    /// Parse a number as f64
    fn parse_number(&mut self, tok: &str) -> Option<Token<'a>> {
        let mut has_dot = false;
        let mut len = 1; // starts with '0'..='9' which has length 1

        for ch in tok.chars().skip(1) {
            if ch == '.' && !has_dot {
                has_dot = true;
            } else if !matches!(ch, '0'..='9') {
                break;
            }
            self.chars.next();
            len += 1; // only accepting '0'..='9' | '.', so char len always == 1
        }

        let tok = &tok[0..len];
        let number = match tok.parse() {
            Ok(a) => a,
            Err(e) => {
                eprintln!("Internal error parsing number: {:?}: {}", tok, e);
                return None;
            }
        };

        Some(Token::Number(number))
    }

    /// Parse an identifier
    fn parse_identifier(&mut self, tok: &'a str, ch: char) -> Option<Token<'a>> {
        let mut len = ch.len_utf8();

        for ch in tok.chars().skip(1) {
            if !ch.is_alphanumeric() {
                break;
            }

            len += ch.len_utf8();
            self.chars.next();
        }

        let ident = &tok[0..len];

        Some(match ident {
            "def" => Token::Def,
            "extern" => Token::Extern,
            _ => Token::Identifier(ident),
        })
    }
}
