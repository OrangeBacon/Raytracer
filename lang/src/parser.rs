use std::borrow::Cow;

use crate::{
    ast::{Error, Expression, Prototype, TopLevel},
    scanner::{NameId, Token, TokenVisitor},
    Scanner,
};

/// Parse a string into an ast
pub struct Parser<'a> {
    /// The source
    scanner: Scanner<'a>,

    /// Has an error token been created?
    has_error: bool,

    /// Is the parser currently recovering from an error
    panic_mode: bool,
}

/// The precedence order for operators
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
enum Precedence {
    #[default]
    None,
    Comparison,
    Multiply,
    Addition,
    Primary,
}

impl<'a> Parser<'a> {
    /// Create a new parser
    pub fn new(scanner: Scanner<'a>) -> Self {
        Self {
            scanner,
            has_error: false,
            panic_mode: false,
        }
    }

    /// Get the next top level ast from the source
    pub fn next(&mut self) -> Option<TopLevel> {
        loop {
            self.synchronise();

            break Some(match self.scanner.peek()? {
                Token::Def => self.def(),
                Token::Extern => self.extern_(),
                Token::Char(';') => continue, // skip top level semicolons
                Token::Identifier(_) | Token::Number(_) | Token::Char(_) => {
                    TopLevel::Expression(self.precedence(Precedence::Comparison))
                }
            });
        }
    }

    /// Synchronise the parser so it can attempt to parse a new statement after
    /// it encounters an error
    fn synchronise(&mut self) {
        if !self.panic_mode {
            return;
        }
        self.panic_mode = false;

        while let Some(next) = self.scanner.peek() {
            match next {
                Token::Def | Token::Extern => return,
                Token::Char(';') => {
                    self.scanner.next();
                    return;
                }
                _ => (),
            }
        }
    }

    /// Parse a function declaration
    fn def(&mut self) -> TopLevel {
        self.scanner.next(); // skip `def` token
        let proto = self.prototype();
        let body = self.precedence(Precedence::Comparison);
        TopLevel::FunctionDefinition(proto, body)
    }

    fn extern_(&mut self) -> TopLevel {
        self.scanner.next(); // skip `extern` token
        let proto = self.prototype();
        TopLevel::ExternDecl(proto)
    }

    /// Parse a function prototype
    fn prototype(&mut self) -> Prototype {
        let name = match self.scanner.next() {
            Some(Token::Identifier(name)) => name,
            Some(tok) => return self.error_ctx("Expected function name", &[tok]),
            None => return self.error("Expected function name, found end of file"),
        };

        if let Some(err) = self.consume(Token::Char('(')) {
            return err;
        }

        let mut args = vec![];
        while let Some(next) = self.scanner.peek() {
            if let Token::Identifier(ident) = next {
                args.push(ident);
                self.scanner.next();
            } else {
                break;
            }

            match self.scanner.peek() {
                Some(Token::Char(',')) => (),
                _ => break,
            }
        }

        if let Some(err) = self.consume(Token::Char(')')) {
            return err;
        }

        Prototype::Prototype { name, args }
    }

    /// Parse an expression with the given precedence
    fn precedence(&mut self, precedence: Precedence) -> Expression {
        let tok = match self.scanner.next() {
            Some(tok) => tok,
            None => return self.error("Expected expression, found end of file"),
        };

        let prefix = tok.visit(Prefix(self));

        let mut expr = match prefix {
            Ok(pre) => pre,
            Err(_) => return self.error_ctx("Expected expression", &[tok]),
        };

        while let Some(peek) = self.scanner.peek() {
            if precedence > peek.visit(TokPrecedence(self)).unwrap_or_default() {
                break;
            }

            let tok = self.scanner.next().unwrap();

            let infix = tok.visit(Infix {
                parser: self,
                lhs: expr,
            });

            match infix {
                Ok(e) => expr = e,
                Err(infix) => {
                    expr = infix.lhs;
                    break;
                }
            };
        }

        expr
    }

    /// Consume a token, otherwise produce an error
    fn consume<T: Error>(&mut self, tok: Token) -> Option<T> {
        Some(match self.scanner.next() {
            Some(t) if t == tok => return None,
            Some(tok) => {
                self.error_ctx(format!("Expected `{}`", self.scanner.fmt_tok(tok)), &[tok])
            }
            None => self.error(format!("Expected open parentheses, `(`, found end of file")),
        })
    }

    /// Trigger an error message while parsing
    fn error<T: Error>(&mut self, message: impl Into<Cow<'static, str>>) -> T {
        self.panic_mode = true;
        self.has_error = true;
        T::new(message)
    }

    /// Trigger an error message while parsing
    fn error_ctx<T: Error>(&mut self, message: impl Into<Cow<'static, str>>, ctx: &[Token]) -> T {
        self.panic_mode = true;
        self.has_error = true;
        T::with_context(message, ctx)
    }
}

struct Prefix<'a, 'b>(&'b mut Parser<'a>);
impl<'a, 'b> TokenVisitor<Expression> for Prefix<'a, 'b> {
    fn identifier(self, name: NameId) -> Result<Expression, Prefix<'a, 'b>> {
        Ok(Expression::Variable(name))
    }

    fn number(self, num: f64) -> Result<Expression, Prefix<'a, 'b>> {
        Ok(Expression::Number(num))
    }
}

struct Infix<'a, 'b> {
    parser: &'b mut Parser<'a>,
    lhs: Expression,
}
impl<'a, 'b> TokenVisitor<Expression> for Infix<'a, 'b> {
    fn char(self, char: char) -> Result<Expression, Infix<'a, 'b>> {
        let prec = Token::Char(char)
            .visit(TokPrecedence(self.parser))
            .unwrap_or_default();
        let prec = prec.next();

        match char {
            '+' => todo!(),
            '-' => todo!(),
            '*' => todo!(),
            '/' => todo!(),
            '<' => todo!(),
            '>' => todo!(),
            _ => Err(self),
        }
    }
}

struct TokPrecedence<'a, 'b>(&'b mut Parser<'a>);
impl<'a, 'b> TokenVisitor<Precedence> for TokPrecedence<'a, 'b> {
    fn char(self, char: char) -> Result<Precedence, TokPrecedence<'a, 'b>> {
        Ok(match char {
            '+' | '-' => Precedence::Addition,
            '*' | '/' => Precedence::Multiply,
            '<' | '>' => Precedence::Comparison,
            _ => Precedence::None,
        })
    }
}

impl Precedence {
    fn next(self) -> Self {
        use Precedence::*;

        match self {
            None => Comparison,
            Comparison => Multiply,
            Multiply => Addition,
            Addition => Primary,
            Primary => Primary,
        }
    }
}
