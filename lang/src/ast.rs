use std::borrow::Cow;

use crate::scanner::{NameId, Token};

/// All top level items in a source file
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum TopLevel {
    FunctionDefinition(Prototype, Expression),
    ExternDecl(Prototype),
    Expression(Expression),
}

/// A function name and arguments
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum Prototype {
    Prototype {
        name: NameId,
        args: Vec<NameId>,
    },
    Error {
        message: Cow<'static, str>,
        tokens: Vec<Token>,
    },
}

/// A single expression tree
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum Expression {
    Number(f64),
    Variable(NameId),
    Binary {
        left: Box<Expression>,
        op: Operator,
        right: Box<Expression>,
    },
    Call {
        callee: Box<Expression>,
        args: Vec<Expression>,
    },
    Error {
        message: Cow<'static, str>,
        tokens: Vec<Token>,
    },
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum Operator {
    Add,
    Subtract,
    Multiply,
    Divide,
    LessThan,
    GreaterThan,
    Error(char),
}

/// Any part of the AST that can turn into an error message
pub trait Error {
    /// Construct a new error message
    fn new(message: impl Into<Cow<'static, str>>) -> Self;

    /// Construct a new error message with some tokens of context
    fn with_context(message: impl Into<Cow<'static, str>>, ctx: &[Token]) -> Self;
}

impl Error for Prototype {
    fn new(message: impl Into<Cow<'static, str>>) -> Self {
        Prototype::Error {
            message: message.into(),
            tokens: vec![],
        }
    }

    fn with_context(message: impl Into<Cow<'static, str>>, ctx: &[Token]) -> Self {
        Prototype::Error {
            message: message.into(),
            tokens: ctx.into(),
        }
    }
}

impl Error for Expression {
    fn new(message: impl Into<Cow<'static, str>>) -> Self {
        Expression::Error {
            message: message.into(),
            tokens: vec![],
        }
    }

    fn with_context(message: impl Into<Cow<'static, str>>, ctx: &[Token]) -> Self {
        Expression::Error {
            message: message.into(),
            tokens: ctx.into(),
        }
    }
}
