use std::borrow::Cow;

use crate::{
    ast::{Error, Expression, Operator, Prototype, TopLevel},
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
    Assignment,
    Or,
    And,
    Equality,
    Comparison,
    Addition,
    Multiply,
    Unary,
    Call,
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
                Token::Char(';') => {
                    self.scanner.next(); // skip top level semicolons
                    continue;
                }
                Token::Identifier(_) | Token::Number(_) | Token::Char(_) => {
                    TopLevel::Expression(self.precedence(Precedence::Comparison))
                }
            });
        }
    }

    /// Print all ast elements to stdout
    pub fn print_ast(mut self) {
        while let Some(node) = self.next() {
            println!("{node:#?}");
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
        let body = self.precedence(Precedence::None);
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
            if precedence > peek.visit(TokPrecedence).unwrap_or_default() {
                break;
            }

            let infix = peek.visit(Infix {
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

    fn char(self, char: char) -> Result<Expression, Self> {
        if char == '(' {
            Ok(self.0.precedence(Precedence::Comparison))
        } else {
            Err(self)
        }
    }
}

struct Infix<'a, 'b> {
    parser: &'b mut Parser<'a>,
    lhs: Expression,
}
impl<'a, 'b> TokenVisitor<Expression> for Infix<'a, 'b> {
    fn char(self, char: char) -> Result<Expression, Infix<'a, 'b>> {
        if char == '(' {
            self.parser.scanner.next(); // consume the operator

            // parse function call
            let mut args = vec![];
            while let Some(_) = self.parser.scanner.peek() {
                args.push(self.parser.precedence(Precedence::Comparison));
                if self.parser.scanner.peek() == Some(Token::Char(')')) {
                    break;
                }
            }

            match self.parser.consume(Token::Char(')')) {
                Some(err) => args.push(err),
                None => (),
            }

            return Ok(Expression::Call {
                callee: Box::new(self.lhs),
                args,
            });
        }

        let operator = match char {
            '+' => Operator::Add,
            '-' => Operator::Subtract,
            '*' => Operator::Multiply,
            '/' => Operator::Divide,
            '<' => Operator::LessThan,
            '>' => Operator::GreaterThan,
            _ => return Err(self),
        };

        self.parser.scanner.next(); // consume the operator

        let prec = Token::Char(char).visit(TokPrecedence).unwrap_or_default();
        let prec = prec.next();
        let rhs = self.parser.precedence(prec);

        Ok(Expression::Binary {
            left: Box::new(self.lhs),
            op: operator,
            right: Box::new(rhs),
        })
    }
}

struct TokPrecedence;
impl TokenVisitor<Precedence> for TokPrecedence {
    fn char(self, char: char) -> Result<Precedence, TokPrecedence> {
        Ok(match char {
            '+' | '-' => Precedence::Addition,
            '*' | '/' => Precedence::Multiply,
            '<' | '>' => Precedence::Comparison,
            '(' => Precedence::Call,
            _ => Precedence::None,
        })
    }
}

impl Precedence {
    fn next(self) -> Self {
        use Precedence::*;

        match self {
            None => Assignment,
            Assignment => Or,
            Or => And,
            And => Equality,
            Equality => Comparison,
            Comparison => Addition,
            Addition => Multiply,
            Multiply => Unary,
            Unary => Call,
            Call => Primary,
            Primary => Primary,
        }
    }
}
