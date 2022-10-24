mod ast;
mod transform_set;

pub use ast::*;

use std::{
    collections::HashSet,
    fmt::Display,
    path::{Path, PathBuf},
    str::{Chars, FromStr},
};

use geometry::{Bounds2, Number, Point2, Point3, Vector3};

use self::transform_set::TransformManager;

/// Error messages encountered while parsing a pbrt file
#[derive(Debug)]
pub enum PbrtError {
    IOError(std::io::Error),
    NonTerminatedString,
    UnexpectedEscapeSequence,
    EofInsideEscapeSequence,
    InvalidNumericConstant,
    InvalidBooleanConstant,
    UnexpectedToken(String),
    UnknownDirectives(Vec<String>),
    EofInsideDirective,
    InvalidKind(&'static str, String),
    InvalidParameter(String, String),
    ExpectedParameter(String),
    UnknownCoordinateSystem(String),
}

impl From<std::io::Error> for PbrtError {
    fn from(err: std::io::Error) -> Self {
        Self::IOError(err)
    }
}

impl<'a, T: Number> From<Token<'a, T>> for PbrtError {
    fn from(tok: Token<'a, T>) -> Self {
        Self::UnexpectedToken(format!("{:?}", tok))
    }
}

impl Display for PbrtError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PbrtError::IOError(err) => {
                write!(f, "IO Error while parsing pbrt file: {err}")
            }
            PbrtError::NonTerminatedString => f.write_str("Non-terminated string"),
            PbrtError::UnexpectedEscapeSequence => f.write_str("Unexpected escape sequence"),
            PbrtError::EofInsideEscapeSequence => {
                f.write_str("Unexpected end of file inside escape sequence")
            }
            PbrtError::InvalidNumericConstant => f.write_str("Unable to parse number"),
            PbrtError::InvalidBooleanConstant => f.write_str("Unable to parse bool"),
            PbrtError::UnexpectedToken(s) => write!(f, "Unexpected token {s:?}"),
            PbrtError::UnknownDirectives(dir) => {
                f.write_str("Unknown directives: ")?;
                let mut set = f.debug_set();
                for dir in dir {
                    set.entry(dir);
                }
                set.finish()
            }
            PbrtError::EofInsideDirective => f.write_str("Unexpected end of file inside directive"),
            PbrtError::InvalidKind(s, k) => write!(f, "{k} is not a recognised kind of {s}"),
            PbrtError::InvalidParameter(ty, name) => {
                write!(f, "parameter {name} of type {ty} unexpected")
            }
            PbrtError::ExpectedParameter(s) => {
                write!(f, "Expected string of \"type name\", got {s:?}")
            }
            PbrtError::UnknownCoordinateSystem(s) => write!(f, "Coordinate system {s:?} unknown"),
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
        let data = std::fs::read_to_string(file_name)?;

        let parser = Parser::new(&data);

        parser.parse()
    }
}

/// State for getting tokens from a string
struct Tokeniser<'a, T: Number> {
    /// The input string char iterator
    chars: Chars<'a>,

    /// Stored data e.g. identifiers, decoded strings, ...
    string: String,

    /// The token that was peeked
    data: Option<TokeniserState<T>>,
}

#[derive(Clone, Copy)]
enum TokeniserState<T: Number> {
    String,
    Identifier,
    Number(T),
    OpenSquareBracket,
    CloseSquareBracket,
}

#[derive(Debug)]
enum Token<'a, T: Number> {
    String(&'a str),
    Number(T),
    Identifier(&'a str),
    OpenSquareBracket,
    CloseSquareBracket,
}

impl<'a, T: Number> Tokeniser<'a, T> {
    /// Create a new tokeniser for the given string
    fn new(data: &'a str) -> Self {
        Self {
            chars: data.chars(),
            string: String::new(),
            data: None,
        }
    }

    /// Try to get the next token in the stream but don't consume it
    fn peek(&mut self) -> Result<Option<Token<T>>, PbrtError> {
        Ok(match self.data {
            Some(data) => Some(self.data_to_token(data)),
            None => {
                self.token()?;
                self.data.map(|data| self.data_to_token(data))
            }
        })
    }

    /// Return and consume the next token in the stream
    fn next(&mut self) -> Result<Option<Token<T>>, PbrtError> {
        Ok(match self.data.take() {
            Some(data) => Some(self.data_to_token(data)),
            None => {
                self.token()?;
                self.data.take().map(|data| self.data_to_token(data))
            }
        })
    }

    /// Convert the internal state into an actual token
    fn data_to_token(&'a self, data: TokeniserState<T>) -> Token<'a, T> {
        match data {
            TokeniserState::String => Token::String(&self.string),
            TokeniserState::Identifier => Token::Identifier(&self.string),
            TokeniserState::Number(n) => Token::Number(n),
            TokeniserState::OpenSquareBracket => Token::OpenSquareBracket,
            TokeniserState::CloseSquareBracket => Token::CloseSquareBracket,
        }
    }

    /// replace the stored token data with new data
    fn token(&mut self) -> Result<(), PbrtError> {
        loop {
            let next = match self.chars.next() {
                Some(c) => c,
                None => return Ok(()),
            };

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
                    self.string.clear();
                    for ch in &mut self.chars {
                        if ch == '"' {
                            self.data = Some(TokeniserState::String);
                            return Ok(());
                        }
                        if ch == '\\' {
                            break;
                        }
                        self.string.push(ch);
                    }

                    while let Some(ch) = self.chars.next() {
                        match ch {
                            '\\' => {
                                let ch = match self.chars.next() {
                                    Some('b') => '\x08',
                                    Some('f') => '\x0C',
                                    Some('n') => '\n',
                                    Some('r') => '\r',
                                    Some('t') => '\t',
                                    Some('\\') => '\\',
                                    Some('\'') => '\'',
                                    Some('\"') => '\"',
                                    Some(_) => return Err(PbrtError::UnexpectedEscapeSequence),
                                    None => return Err(PbrtError::EofInsideEscapeSequence),
                                };

                                self.string.push(ch);
                            }
                            '"' => {
                                self.data = Some(TokeniserState::String);
                                return Ok(());
                            }
                            _ => self.string.push(ch),
                        }
                    }
                    return Err(PbrtError::NonTerminatedString);
                }
                '[' => {
                    self.data = Some(TokeniserState::OpenSquareBracket);
                    return Ok(());
                }
                ']' => {
                    self.data = Some(TokeniserState::CloseSquareBracket);
                    return Ok(());
                }
                _ => {
                    self.string.clear();
                    self.string.push(next);
                    for ch in self.chars.clone() {
                        if " \n\t\r\"[]".contains(ch) {
                            break;
                        }
                        self.chars.next();
                        self.string.push(ch);
                    }

                    if "0123456789.-".contains(next) {
                        let num = self
                            .string
                            .parse()
                            .map_err(|_| PbrtError::InvalidNumericConstant)?;
                        self.data = Some(TokeniserState::Number(num));
                        return Ok(());
                    }

                    self.data = Some(TokeniserState::Identifier);
                    return Ok(());
                }
            }
        }
    }
}

macro_rules! expect {
    ($tokeniser:expr, $matcher:ident) => {
        match $tokeniser.next()? {
            Some(Token::$matcher) => (),
            Some(tok) => return Err(PbrtError::UnexpectedToken(format!("{:?}", tok))),
            None => return Err(PbrtError::EofInsideDirective),
        }
    };
    ($tokeniser:expr, $matcher:ident()) => {
        match $tokeniser.next()? {
            Some(Token::$matcher(val)) => val,
            Some(tok) => return Err(PbrtError::UnexpectedToken(format!("{:?}", tok))),
            None => return Err(PbrtError::EofInsideDirective),
        }
    };
}

/// Macro to make it easier to match directive parameters to their fields in ast
macro_rules! directive {
    (match $tokeniser:expr => $(
        $pattern:pat => $(@$modifier:ident)? $path:expr
    ),* $(,)?) => {
        while let Some(peek) = $tokeniser.peek()? {
            // break if there is an identifier because it probably starts a new directive
            match peek {
                Token::String(_) => (),
                Token::Identifier(_) => break,
                _ => return Err(PbrtError::UnexpectedToken(format!("{peek:?}"))),
            }

            // this should always succeed due to checking the peeked value
            let ident = match $tokeniser.next()? {
                Some(Token::String(s)) => s,
                Some(tok) => return Err(PbrtError::UnexpectedToken(format!("{tok:?}"))),
                None => return Err(PbrtError::EofInsideDirective),
            };
            match ident.split_once(' ') {
                $(
                    Some($pattern) => directive!(@set $tokeniser => $(@$modifier)? $path),
                )*

                Some((ty, name)) => {
                    return Err(PbrtError::InvalidParameter(
                        ty.to_string(),
                        name.to_string(),
                    ))
                }
                None => return Err(PbrtError::ExpectedParameter(ident.to_string())),
            }
        }
    };
    (@set $tokeniser:expr => $path:expr) => {
        $path = <_>::parse(&mut $tokeniser)?
    };
    (@set $tokeniser:expr => @num $path:expr) => {
        $path = <NumberWrapper<_>>::parse(&mut $tokeniser)?.0
    };
    (@set $tokeniser:expr => @optnum $path:expr) => {
        $path = <Option<NumberWrapper<_>>>::parse(&mut $tokeniser)?.map(|d| d.0)
    };
}

struct Parser<'a, T: Number> {
    tokeniser: Tokeniser<'a, T>,
    result_file: PbrtFile<T>,
    transform: TransformManager<T>,
}

impl<'a, T: Number> Parser<'a, T> {
    fn new(data: &'a str) -> Self {
        Self {
            tokeniser: Tokeniser::new(data),
            result_file: PbrtFile::default(),
            transform: TransformManager::default(),
        }
    }

    fn parse(mut self) -> Result<PbrtFile<T>, PbrtError> {
        let mut unknown_directives = HashSet::new();
        let mut skip_unknown = false;

        while let Some(tok) = self.tokeniser.next()? {
            match tok {
                Token::Identifier("WorldBegin") => break, // only looking at global options atm
                Token::Identifier(ident) => {
                    skip_unknown = false;
                    let directive = ident.parse();
                    let directive = match directive {
                        Ok(dir) => dir,
                        Err(PbrtError::UnknownDirectives(dir)) => {
                            for dir in dir {
                                unknown_directives.insert(dir);
                            }
                            skip_unknown = true;
                            continue;
                        }
                        Err(err) => return Err(err),
                    };
                    match self.parse_directive(directive) {
                        Ok(()) => (),
                        Err(PbrtError::UnknownDirectives(dir)) => {
                            for dir in dir {
                                unknown_directives.insert(dir);
                            }
                            skip_unknown = true;
                            continue;
                        }
                        Err(err) => return Err(err),
                    };
                }
                tok => {
                    if skip_unknown {
                        continue;
                    }
                    return Err(tok.into());
                }
            }
        }

        if !unknown_directives.is_empty() {
            println!("Warning: Skipping unknown directives: {unknown_directives:?}",);
        }

        Ok(self.result_file)
    }

    fn parse_directive(&mut self, directive: Directives) -> Result<(), PbrtError> {
        match directive {
            // TRANSFORMS
            Directives::Identity => self.transform.identity(),
            Directives::Translate => {
                let [dx, dy, dz] = self.float_list()?;
                self.transform.translate(dx, dy, dz);
            }
            Directives::Scale => {
                let [sx, sy, sz] = self.float_list()?;
                self.transform.scale(sx, sy, sz);
            }
            Directives::Rotate => {
                let [angle, ax, ay, az] = self.float_list()?;
                self.transform.rotate(angle, ax, ay, az);
            }
            Directives::LookAt => {
                let pos = Point3::from_array(self.float_list()?);
                let look = Point3::from_array(self.float_list()?);
                let up = Vector3::from_array(self.float_list()?);
                self.transform.look_at(pos, look, up);
            }
            Directives::CoordinateSystem => {
                let name = expect!(self.tokeniser, String());
                self.transform.insert(name);
            }
            Directives::CoordSysTransform => {
                let name = expect!(self.tokeniser, String());
                let trans = match self.transform.get(name) {
                    Some(t) => t,
                    None => return Err(PbrtError::UnknownCoordinateSystem(name.to_string())),
                };
                self.transform.set(trans);
            }
            Directives::Transform => {
                let transform = self.float_list()?;
                self.transform.set_transform(transform);
            }
            Directives::ConcatTransform => {
                let transform = self.float_list()?;
                self.transform.concat(transform);
            }

            // GLOBAL SETTINGS
            Directives::Camera => match expect!(self.tokeniser, String()) {
                "environment" => {
                    let mut cam = CameraEnvironment::default();
                    directive! { match self.tokeniser =>
                        ("float", "shutteropen") => @num self.result_file.camera.shutter_open,
                        ("float", "shutterclose") => @num self.result_file.camera.shutter_close,
                        ("float", "screenwindow") => cam.screen_window,
                        ("float", "frameaspectratio") => @optnum cam.frame_aspect_ratio,
                    }
                    self.result_file.camera.kind = CameraKind::Environment(cam);
                    self.result_file.camera.transform = self.transform.insert("camera");
                }
                "orthographic" => {
                    let mut cam = CameraOrthographic::default();
                    directive! { match self.tokeniser =>
                        ("float", "shutteropen") => @num self.result_file.camera.shutter_open,
                        ("float", "shutterclose") => @num self.result_file.camera.shutter_close,
                        ("float", "screenwindow") => cam.screen_window,
                        ("float", "frameaspectratio") => @optnum cam.frame_aspect_ratio,
                        ("float", "lensradius") => @num cam.lens_radius,
                        ("float", "focaldistance") => @num cam.focal_distance,
                    }
                    self.result_file.camera.kind = CameraKind::Orthographic(cam);
                    self.result_file.camera.transform = self.transform.insert("camera");
                }
                "perspective" => {
                    let mut cam = CameraPerspective::default();
                    directive! { match self.tokeniser =>
                        ("float", "shutteropen") => @num self.result_file.camera.shutter_open,
                        ("float", "shutterclose") => @num self.result_file.camera.shutter_close,
                        ("float", "screenwindow") => cam.screen_window,
                        ("float", "frameaspectratio") => @optnum cam.frame_aspect_ratio,
                        ("float", "lensradius") => @num cam.lens_radius,
                        ("float", "focaldistance") => @num cam.focal_distance,
                        ("float", "fov") => @num cam.fov,
                        ("float", "halffov") => @optnum cam.half_fov,
                    }
                    self.result_file.camera.kind = CameraKind::Perspective(cam);
                    self.result_file.camera.transform = self.transform.insert("camera");
                }
                "realistic" => {
                    let mut cam = CameraRealistic::default();
                    directive! { match self.tokeniser =>
                        ("float", "shutteropen") => @num self.result_file.camera.shutter_open,
                        ("float", "shutterclose") => @num self.result_file.camera.shutter_close,
                        ("string", "lensfile") => cam.lens_file,
                        ("float", "aperturediameter") => @num cam.aperture_diameter,
                        ("float", "focusdistance") => @num cam.focus_distance,
                        ("bool", "simpleweighting") => cam.simple_weighting
                    }
                    self.result_file.camera.kind = CameraKind::Realistic(cam);
                    self.result_file.camera.transform = self.transform.insert("camera");
                }
                str => return Err(PbrtError::InvalidKind("Camera", str.to_string())),
            },
            Directives::Sampler => match expect!(self.tokeniser, String()) {
                "stratified" => {
                    let mut sampler = SamplerStratified::default();
                    directive! { match self.tokeniser =>
                        ("bool", "jitter") => sampler.jitter,
                        ("integer", "xsamples") => sampler.x_samples,
                        ("integer", "ysamples") => sampler.y_samples,
                    }
                    self.result_file.sampler = Sampler::Stratified(sampler);
                }
                str => {
                    if let Ok(kind) = str.parse() {
                        let mut sampler = SamplerPixel {
                            kind,
                            ..Default::default()
                        };
                        directive! { match self.tokeniser =>
                            ("integer", "pixelsamples") => sampler.pixel_samples,
                        }
                        self.result_file.sampler = Sampler::Pixel(sampler);
                    } else {
                        return Err(PbrtError::InvalidKind("Sampler", str.to_string()));
                    }
                }
            },
            Directives::Film => match expect!(self.tokeniser, String()) {
                "image" => {
                    self.result_file.film = Default::default();
                    directive! { match self.tokeniser =>
                        ("integer", "xresolution") => self.result_file.film.x_resolution,
                        ("integer", "yresolution") => self.result_file.film.y_resolution,
                        ("float", "cropwindow") => self.result_file.film.crop_window,
                        ("float", "scale") => @num self.result_file.film.scale,
                        ("float", "maxsampleluminance") => @num self.result_file.film.max_sample_luminance,
                        ("float", "diagonal") => @num self.result_file.film.diagonal,
                        ("string", "filename") => self.result_file.film.file_name,
                    }
                }
                str => return Err(PbrtError::InvalidKind("Film", str.to_string())),
            },
            Directives::PixelFilter => match expect!(self.tokeniser, String()) {
                "box" => {
                    let mut filter = FilterBox::default();
                    directive! { match self.tokeniser =>
                        ("float", "xwidth") => @num filter.x_width,
                        ("float", "ywidth") => @num filter.y_width,
                    }
                    self.result_file.filter = Filter::Box(filter);
                }
                "gaussian" => {
                    let mut filter = FilterGaussian::default();
                    directive! { match self.tokeniser =>
                        ("float", "xwidth") => @num filter.x_width,
                        ("float", "ywidth") => @num filter.y_width,
                        ("float", "alpha") => @num filter.alpha,
                    }
                    self.result_file.filter = Filter::Gaussian(filter);
                }
                "mitchell" => {
                    let mut filter = FilterMitchell::default();
                    directive! { match self.tokeniser =>
                        ("float", "xwidth") => @num filter.x_width,
                        ("float", "ywidth") => @num filter.y_width,
                        ("float", "B") => @num filter.b,
                        ("float", "C") => @num filter.c,
                    }
                    self.result_file.filter = Filter::Mitchell(filter);
                }
                "sinc" => {
                    let mut filter = FilterSinc::default();
                    directive! { match self.tokeniser =>
                        ("float", "xwidth") => @num filter.x_width,
                        ("float", "ywidth") => @num filter.y_width,
                        ("float", "tau") => @num filter.tau,
                    }
                    self.result_file.filter = Filter::Sinc(filter);
                }
                "triangle" => {
                    let mut filter = FilterTriangle::default();
                    directive! { match self.tokeniser =>
                        ("float", "xwidth") => @num filter.x_width,
                        ("float", "ywidth") => @num filter.y_width,
                    }
                    self.result_file.filter = Filter::Triangle(filter);
                }
                str => return Err(PbrtError::InvalidKind("Filter", str.to_string())),
            },

            tok => return Err(PbrtError::UnknownDirectives(vec![format!("{tok:?}")])),
        }
        Ok(())
    }

    /// Get exactly N floats
    fn float_list<const N: usize>(&mut self) -> Result<[T; N], PbrtError> {
        let mut result = [T::ZERO; N];

        for elem in &mut result {
            *elem = expect!(self.tokeniser, Number());
        }

        Ok(result)
    }
}

/// All supported directives
#[derive(Debug)]
enum Directives {
    // transforms
    Identity,
    Translate,
    Scale,
    Rotate,
    LookAt,
    CoordinateSystem,
    CoordSysTransform,
    Transform,
    ConcatTransform,

    // setup
    Camera,
    Sampler,
    Film,
    PixelFilter,
    Accelerator,
    Integrator,
    TransformTimes,

    // both
    MakeNamedMedium, // warn animated
    MediumInterface,

    // world description
    AttributeBegin,
    AttributeEnd,
    ActiveTransform,
    AreaLightSource,
    LightSource,       // warn animated
    MakeNamedMaterial, // warn animated
    NamedMaterial,
    ObjectBegin,
    ObjectInstance,
    ReverseOrientation,
    Shape,
    TransformBegin,
    TransformEnd,
    Texture, // warn animated
}

impl FromStr for Directives {
    type Err = PbrtError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use Directives::*;
        Ok(match s {
            "Identity" => Identity,
            "Translate" => Translate,
            "Scale" => Scale,
            "Rotate" => Rotate,
            "LookAt" => LookAt,
            "CoordinateSystem" => CoordinateSystem,
            "CoordSysTransform" => CoordSysTransform,
            "Transform" => Transform,
            "ConcatTransform" => ConcatTransform,
            "Camera" => Camera,
            "Film" => Film,
            "Sampler" => Sampler,
            "PixelFilter" => PixelFilter,
            "Accelerator" => Accelerator,
            "Integrator" => Integrator,
            "TransformTimes" => TransformTimes,
            "MakeNamedMedium" => MakeNamedMedium,
            "MediumInterface" => MediumInterface,
            "AttributeBegin" => AttributeBegin,
            "AttributeEnd" => AttributeEnd,
            "ActiveTransform" => ActiveTransform,
            "AreaLightSource" => AreaLightSource,
            "LightSource" => LightSource,
            "MakeNamedMaterial" => MakeNamedMaterial,
            "NamedMaterial" => NamedMaterial,
            "ObjectBegin" => ObjectBegin,
            "ObjectInstance" => ObjectInstance,
            "ReverseOrientation" => ReverseOrientation,
            "Shape" => Shape,
            "TransformBegin" => TransformBegin,
            "TransformEnd" => TransformEnd,
            "Texture" => Texture,
            _ => return Err(PbrtError::UnknownDirectives(vec![s.into()])),
        })
    }
}

impl FromStr for SamplerPixelKind {
    type Err = PbrtError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use SamplerPixelKind::*;
        Ok(match s {
            "02sequence" | "lowdiscrepancy" => ZeroTwoSequence,
            "halton" => Halton,
            "maxmindist" => MaxMinDist,
            "random" => Random,
            "sobol" => Sobol,
            _ => return Err(PbrtError::InvalidKind("Sampler", s.to_string())),
        })
    }
}

/// Any value that can be parsed as a parameter value from a token stream
trait ParseValue<T: Number>
where
    Self: Sized,
{
    fn parse(tokeniser: &mut Tokeniser<T>) -> Result<Self, PbrtError>;
}

impl<T: Number, U: ParseValue<T>> ParseValue<T> for Option<U> {
    fn parse(tokeniser: &mut Tokeniser<T>) -> Result<Self, PbrtError> {
        Ok(Some(U::parse(tokeniser)?))
    }
}

/// Wrapper required due to orphan rules not letting there be a blanket impl for
/// `T: Number`.  This therefore makes it require the special cases in the
/// `directive!` macro.
struct NumberWrapper<T: Number>(T);
impl<T: Number> ParseValue<T> for NumberWrapper<T> {
    fn parse(tokeniser: &mut Tokeniser<T>) -> Result<Self, PbrtError> {
        let must_close = match tokeniser.peek()? {
            Some(Token::OpenSquareBracket) => {
                tokeniser.next()?;
                true
            }
            _ => false,
        };

        let num = expect!(tokeniser, Number());

        if must_close {
            expect!(tokeniser, CloseSquareBracket);
        }

        Ok(NumberWrapper(num))
    }
}

impl<T: Number> ParseValue<T> for Bounds2<T> {
    fn parse(tokeniser: &mut Tokeniser<T>) -> Result<Self, PbrtError> {
        expect!(tokeniser, OpenSquareBracket);

        let mut vals = [T::ZERO; 4];
        for val in &mut vals {
            *val = expect!(tokeniser, Number());
        }

        expect!(tokeniser, CloseSquareBracket);

        Ok(Bounds2::new(
            Point2::new(vals[0], vals[1]),
            Point2::new(vals[2], vals[3]),
        ))
    }
}

impl<T: Number> ParseValue<T> for PathBuf {
    fn parse(tokeniser: &mut Tokeniser<T>) -> Result<Self, PbrtError> {
        let path = expect!(tokeniser, String());

        Ok(path.into())
    }
}

impl<T: Number> ParseValue<T> for bool {
    fn parse(tokeniser: &mut Tokeniser<T>) -> Result<Self, PbrtError> {
        let bool = expect!(tokeniser, String());

        match bool {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(PbrtError::InvalidBooleanConstant),
        }
    }
}

impl<T: Number> ParseValue<T> for u32 {
    fn parse(tokeniser: &mut Tokeniser<T>) -> Result<Self, PbrtError> {
        let must_close = match tokeniser.peek()? {
            Some(Token::OpenSquareBracket) => {
                tokeniser.next()?;
                true
            }
            _ => false,
        };
        let num = expect!(tokeniser, Number());

        if num % T::ONE != T::ZERO {
            return Err(PbrtError::InvalidNumericConstant);
        }

        if must_close {
            expect!(tokeniser, CloseSquareBracket);
        }

        Ok(num.i32() as _)
    }
}
