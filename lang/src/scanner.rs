use std::{collections::HashMap, fmt::Display, str::Chars};

/// A single token
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum Token {
    Def,
    Extern,
    Identifier(NameId),
    Number(f64),
    Char(char),
}

#[allow(unused_variables)]
pub trait TokenVisitor<T>
where
    Self: Sized,
{
    fn def(self) -> Result<T, Self> {
        Err(self)
    }

    fn extern_(self) -> Result<T, Self> {
        Err(self)
    }

    fn identifier(self, name: NameId) -> Result<T, Self> {
        Err(self)
    }

    fn number(self, num: f64) -> Result<T, Self> {
        Err(self)
    }

    fn char(self, char: char) -> Result<T, Self> {
        Err(self)
    }
}

/// A 'not iterator' to create tokens from a string
#[derive(Debug, Clone)]
pub struct Scanner<'a> {
    chars: Chars<'a>,
    peeked: Option<Token>,
    name_to_id: HashMap<&'a str, NameId>,
    id_to_name: HashMap<NameId, &'a str>,
}

/// Index into the scanner's name table
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct NameId(usize);

impl<'a> Scanner<'a> {
    /// Create a new scanner from the provided source code
    pub fn new(source: &'a str) -> Self {
        Self {
            chars: source.chars(),
            peeked: None,
            name_to_id: HashMap::new(),
            id_to_name: HashMap::new(),
        }
    }

    /// Print all the tokens in the source string to stdout
    pub fn print_tokens(mut self) {
        while let Some(tok) = self.next() {
            println!("{:?}", tok);
        }
    }

    /// Get the next token from the source code
    pub fn next(&mut self) -> Option<Token> {
        match self.peeked.take() {
            Some(tok) => Some(tok),
            None => self.token(),
        }
    }

    /// Peek the next token from the source code without consuming it
    pub fn peek(&mut self) -> Option<Token> {
        match self.peeked {
            Some(tok) => Some(tok),
            None => {
                let tok = self.token();
                self.peeked = tok;
                tok
            }
        }
    }

    /// Parse a new token
    fn token(&mut self) -> Option<Token> {
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
    fn parse_number(&mut self, tok: &str) -> Option<Token> {
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
    fn parse_identifier(&mut self, tok: &'a str, ch: char) -> Option<Token> {
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
            _ => {
                if let Some(&id) = self.name_to_id.get(ident) {
                    Token::Identifier(id)
                } else {
                    let id = NameId(self.name_to_id.len());
                    self.name_to_id.insert(ident, id);
                    self.id_to_name.insert(id, ident);
                    Token::Identifier(id)
                }
            }
        })
    }

    /// Get a formatter for a single token to make displaying it look better
    pub fn fmt_tok(&'a self, tok: Token) -> impl 'a + Display {
        TokenDisplay(self, tok)
    }
}

impl Token {
    /// Call a visitor on a token
    pub fn visit<T, V: TokenVisitor<T>>(&self, v: V) -> Result<T, V> {
        match self {
            Token::Def => v.def(),
            Token::Extern => v.extern_(),
            Token::Identifier(name) => v.identifier(*name),
            Token::Number(num) => v.number(*num),
            Token::Char(char) => v.char(*char),
        }
    }
}

struct TokenDisplay<'a, 'b>(&'a Scanner<'b>, Token);
impl Display for TokenDisplay<'_, '_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.1 {
            Token::Def => f.write_str("def"),
            Token::Extern => f.write_str("extern"),
            Token::Identifier(ident) => write!(f, "Identifier \"{}\"", self.0.id_to_name[&ident]),
            Token::Number(num) => write!(f, "Number {}", num),
            Token::Char(ch) => write!(f, "Character \"{:?}\"", ch),
        }
    }
}
