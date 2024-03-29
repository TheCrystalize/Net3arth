const StandardLibNamesArray =
`identity
draw
blurSphere
blurCube
blurSphereVolume
blurCubeVolume
blurTetrahedron
blurNormalCube
blurOctahedron
blurIcosahedron
blurDodecahedron
blurTeapot
torus
umbilicTorus
stl
mobius3D
hypershift3D
bubble3D
julian3D
juliaq3D
parabola3D
sphereInv
trigCosh3D
trigExp3D
trigLog3D
trigSinh3D
trigTanh3D
unbubble3D
matrix3D
scale3D
scale3D3
translate3D
rotate3D
perspective3D
viewBox
viewSphere
blurImage
setImage
reset
arcsinh
arctanh
bent
blurCircle
blurGasket
blurGaussian
blurNgon
blurSine
blurSquare
blurTriangle
bTransform
bubble
circleInv
cpow
cylinder
dc_poincareDisc
disc
dragon
ePush
eRotate
flipX
flipY
hypershape
hypershift
hypertile3
jac_cd
jac_cn
jac_dn
jac_elk
jac_sn
julian
juliaq
juliascope
mobius
multiMobius
murl2
nSplit
pdj
pointSymmetry
rotate
scale
scale2
schwarzChristoffelMap
schwarzChristoffelInverseMap
schwarzTriangle
sinusoidal
skew
smartcrop
smartshape
splits
tileHelp
tileLog
translate
trigCosh
trigExp
trigLog
trigSinh
trigTanh
unbubble
weierstrassElliptic
brighten
color
dither
shaderPass
normalizeColors
rainbowCirc
rainbowCircAdd
paletteMod
gamma
gradient
repeatingGradient
hslShift
lerpColor
lerpHSL
lerpRGB
normalizeColor
setHue
setAlpha
setSaturation
setLightness
normalMap
matcap
heightMap
mist
basicLighting
basicEnvironment
basicEnvironmentOrth
advancedLighting
advancedLightingOrth
ambientOcclusion
ambientOcclusionBig
ambientOcclusion2
edgeDetection
specular
specularOrth`.split('\n');

const StandardLibHelperNamesArray =
`decodeSTL
rgbToHsl
hslToRgb
conj
div
divScalar
getIrrationals
squircleDistance
add
neg
sub
addScalar
mult
multScalar
sqrt
lerp
log
pow
exp
dot
sin
cos
tan
sinh
cosh
tanh
gaussRnd
jacobiAm
jacobiAmA
jacobiAmM
jacobiAmK
addPoly
multiplyPoly
multiplyMatrices
findRoots
zeta
toRadian
applyMatrix
multiplyMatrixAndPoint
matrixMultiply
crossProduct
dotProduct
lerp4
normalOf3Points
normalize3
matrixTranslate
matrixScale
matrixRotateX
matrixRotateY
matrixRotateZ
matrixPerspectiveProjection
triangle3D
getNormal
zBuffer
firstBuffer
lastBuffer
averageBuffer`.split('\n');

const StandardLibConstructorNamesArray =
`C
colorRGB
colorHSL`.split('\n');


const StandardLibConstantsNamesArray =
`DEGREE
IDENTITY_MATRIX`.split('\n');

let StandardLibNames = StandardLibNamesArray.join('|');

let StandardLibHelperNames = StandardLibHelperNamesArray.join('|');

let StandardLibConstructorNames = StandardLibConstructorNamesArray.join('|');

let StandardLibConstantsNames = StandardLibConstantsNamesArray.join('|');

// helful resource:
// https://github.com/ajaxorg/ace/wiki/Creating-or-Extending-an-Edit-Mode#common-tokens

/* ***** BEGIN LICENSE BLOCK *****
 * Distributed under the BSD license:
 *
 * Copyright (c) 2010, Ajax.org B.V.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of Ajax.org B.V. nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL AJAX.ORG B.V. BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ***** END LICENSE BLOCK ***** */
/*
 * Adapted by the Net3arth team.
 */

define('ace/mode/3arthLang',function(require, exports, module) {
"use strict";

var oop = require("../lib/oop");
var DocCommentHighlightRules = require("./doc_comment_highlight_rules").DocCommentHighlightRules;
var TextHighlightRules = require("./text_highlight_rules").TextHighlightRules;

var identifierRe = "[a-zA-Z\\$_\u00a1-\uffff][a-zA-Z\\d\\$_\u00a1-\uffff]*";

var EarthLangHighlightRules = function(options) {
    // see: https://developer.mozilla.org/en/JavaScript/Reference/Global_Objects
    var keywordMapper = this.createKeywordMapper({
        "support.class":
            "Array|Boolean|Date|Function|Iterator|Number|Object|RegExp|String|Proxy|"  + // Constructors
            "Namespace|QName|XML|XMLList|"                                             + // E4X
            "ArrayBuffer|Float32Array|Float64Array|Int16Array|Int32Array|Int8Array|"   +
            "Uint16Array|Uint32Array|Uint8Array|Uint8ClampedArray|"                    +
            "Error|EvalError|InternalError|RangeError|ReferenceError|StopIteration|"   + // Errors
            "SyntaxError|TypeError|URIError|"                                          +
            "decodeURI|decodeURIComponent|encodeURI|encodeURIComponent|eval|isFinite|" + // Non-constructor functions
            "isNaN|parseFloat|parseInt|"                                               +
            "JSON|Math|oldBuffer|newBuffer|"                                           + // Other
            "this|arguments|prototype|window|document"                                 , // Pseudo
        "keyword":
            "const|yield|import|get|set|async|await|" +
            "break|case|catch|continue|default|delete|do|else|finally|for|function|" +
            "if|in|of|instanceof|new|return|switch|throw|try|typeof|let|var|while|with|debugger|" +
            // invalid or reserved
            "__parent__|__count__|escape|unescape|with|__proto__|" +
            "class|enum|extends|super|export|implements|private|public|interface|package|protected|static",
        "storage.type":
            "const|let|var|function",
        "constant.language":
            "null|Infinity|NaN|undefined|"+StandardLibConstantsNames,
        "support.function":
            "alert",
        "constant.language.boolean": "true|false",
        "support.function": StandardLibHelperNames,
    }, "identifier");

    // keywords which can be followed by regular expressions
    var kwBeforeRe = "case|do|else|finally|in|instanceof|return|throw|try|typeof|yield|void";

    var escapedRe = "\\\\(?:x[0-9a-fA-F]{2}|" + // hex
        "u[0-9a-fA-F]{4}|" + // unicode
        "u{[0-9a-fA-F]{1,6}}|" + // es6 unicode
        "[0-2][0-7]{0,2}|" + // oct
        "3[0-7][0-7]?|" + // oct
        "[4-7][0-7]?|" + //oct
        ".)";
    // regexp must not have capturing parentheses. Use (?:) instead.
    // regexps are ordered -> the first match is used

    this.$rules = {
        "qqstringType" : [
            {
                token : "constant.language.escape",
                regex : escapedRe
            }, {
                token : "string",
                regex : "\\\\$",
                consumeLineEnd  : true
            }, {
                token : "string",
                regex : '"|$',
                next  : "earthLangFunction"
            }, {
                defaultToken: "string"
            }
        ],
        "qstringType" : [
            {
                token : "constant.language.escape",
                regex : escapedRe
            }, {
                token : "string",
                regex : "\\\\$",
                consumeLineEnd  : true
            }, {
                token : "string",
                regex : "'|$",
                next  : "earthLangFunction"
            }, {
                defaultToken: "string"
            }
        ],
        "qqstringParam" : [
            {
                token : "constant.language.escape",
                regex : escapedRe
            }, {
                token : "string",
                regex : "\\\\$",
                consumeLineEnd  : true
            }, {
                token : "string",
                regex : '"|$',
                next  : "earthLangParams"
            }, {
                defaultToken: "string"
            }
        ],
        "qstringParam" : [
            {
                token : "constant.language.escape",
                regex : escapedRe
            }, {
                token : "string",
                regex : "\\\\$",
                consumeLineEnd  : true
            }, {
                token : "string",
                regex : "'|$",
                next  : "earthLangParams"
            }, {
                defaultToken: "string"
            }
        ],


        "earthLangTransformPost": [
            comments("earthLangTransformPost"),
            {
                token: "empty",
                regex: ";",
                next: "pop"
            },
            {
                token: "keyword.control",
                regex: "->",
                next: "earthLangTransform"
            }
        ],
        "earthLangParams": [
            comments("earthLangParams"),
            {
                token : "empty",
                regex : "'\\(|\\)"
            },
            {
                token : "string",
                regex : "'(?=.)",
                next  : "qstringParam"
            },
            {
                token : "string",
                regex : '"(?=.)',
                next  : "qqstringParam"
            },
            {
                token: "empty",
                regex: ";",
                next: "pop"
            },
            {
                token: "keyword.control",
                regex: "->",
                next: "earthLangTransform"
            },
            {
                token: "keyword.control",
                regex: "\\*|\\/|\\+|\\-",
            },
            {
                token: "constant.language",
                regex: "true|false"
            },
            {
                token: "constant.language",
                regex: "i\\b",
            },
            {
                token : "constant.numeric", // hexadecimal, octal and binary
                regex : /0(?:[xX][0-9a-fA-F]+|[oO][0-7]+|[bB][01]+)\b/
            }, {
                token : "constant.numeric", // decimal integers and floats
                regex : /(?:\d\d*(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+\b)?/
            },
            {
                token: "support.class",
                regex: "(?:" + StandardLibConstructorNames + ")\\b"
            },
            {
                token: "entity.function.name",
                regex: "\\w+\\b"
            }
        ],
        "earthLangTransform": [
            comments("earthLangTransform"),
            {
                token: "keyword.control",
                regex: "(choose|sum|sumColor|productColor|product)\\b",
                next: "earthLangChoose"
            },
            {
                token: "keyword.control",
                regex: "(xaos)\\b",
                next: "earthLangXaos"
            },
            {
                token: "keyword.control",
                regex: "(switch)\\b",
                next: "earthLangSwitch"
            },
            {
                token: "variable.language",
                regex: "(?:" + StandardLibNames + ")\\b"
            },
            {
                token: "entity.function.name",
                regex: "\\w+\\b"
            },
            {
                token: "empty",
                regex: "\\(",
                next: "earthLangParams"
            }
        ],

        "earthLangChoose": [
            comments("earthLangChoose"),
            {
                token: "empty",
                regex: "{"
            },
            {
                token: "empty",
                regex: "}",
                next: "earthLangTransformPost"
            },
            {
                token: "keyword.control",
                regex: ":",
                push: "earthLangChoose",
                next: "earthLangTransform"
            },{
                token : "support.class", // hexadecimal, octal and binary
                regex : /0(?:[xX][0-9a-fA-F]+|[oO][0-7]+|[bB][01]+)\b/
            }, {
                token : "support.class", // decimal integers and floats
                regex : /(?:\d\d*(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+\b)?/
            }, {
                token : "keyword.control",
                regex : /--|\+\+|\.{3}|===|==|=|!=|!==|<+=?|>+=?|!|&&|\|\||\?:|[!$%&*+\-~\/^]=?/
            },
        ],

        "earthLangSwitch": [
          comments("earthLangSwitch"),
          {
            token: "empty",
            regex: "{",
          },
          {
            token: "keyword.control",
            regex: "(case)\\b",
          },
          {
            token: "keyword.control",
            regex: "(default)\\b",
          },
          {
            token: "empty",
            regex: "\\(",
            push: "earthLangSwitch",
            next: "startjs"
          },
          {
            token: "empty",
            regex: "}",
            next: "earthLangTransformPost"
          },
          {
            token: "keyword.control",
            regex: ":",
            push: "earthLangSwitch",
            next: "earthLangTransform"
          }
        ],

        "earthLangXaos": [
            comments("earthLangXaos"),
            {
                token: "empty",
                regex: "{"
            },
            {
                token: "empty",
                regex: "}",
                next: "earthLangTransformPost"
            },
            {
                token: "keyword.control",
                regex: ":",
                push: "earthLangXaos",
                next: "earthLangXaosIn"
            },{
                token : "support.class", // hexadecimal, octal and binary
                regex : /0(?:[xX][0-9a-fA-F]+|[oO][0-7]+|[bB][01]+)\b/
            }, {
                token : "support.class", // decimal integers and floats
                regex : /(?:\d\d*(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+\b)?/
            }, {
                token : "keyword.control",
                regex : /--|\+\+|\.{3}|===|==|=|!=|!==|<+=?|>+=?|!|&&|\|\||\?:|[!$%&*+\-~\/^]=?/
            },
        ],
        "earthLangXaosIn": [
            comments("earthLangXaosIn"),
            {
                token: "empty",
                regex: "{"
            },
            {
                token: "empty",
                regex: "}",
                next: "earthLangTransformPost"
            },
            {
                token: "keyword.control",
                regex: ":",
                next: "earthLangXaosMap"
            },
            {
                token: "support.class",
                regex: "e|o|_",
            },
        ],
        "earthLangXaosMap": [
            comments("earthLangXaosMap"),
            {
                token: "empty",
                regex: "{"
            },
            {
                token: "empty",
                regex: "}",
                next: "earthLangTransformPost"
            },
            {
                token: "keyword.control",
                regex: ":",
                next: "earthLangTransform"
            },{
                token : "support.class", // hexadecimal, octal and binary
                regex : /0(?:[xX][0-9a-fA-F]+|[oO][0-7]+|[bB][01]+)\b/
            }, {
                token : "support.class", // decimal integers and floats
                regex : /(?:\d\d*(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+\b)?/
            }, {
                token : "keyword.control",
                regex : /--|\+\+|\.{3}|===|==|=|!=|!==|<+=?|>+=?|!|&&|\|\||\?:|[!$%&*+\-~\/^]=?/
            },
        ],

        "earthLangFunctionType": [
            comments("earthLangFunctionType"),
            {
                token: "storage.type",
                regex: "(?:number|complex|bool|string|array|object|function)",
                push: "earthLangFunctionType",
                next: "earthLangFunction"
            },
            {
                token: "storage.type",
                regex: "\\w+\\b"
            },
            {
                token: "empty",
                regex: "\\(",
            },
            {
                token: "empty",
                regex: "\\)",
                next: "startjs"
            },
        ],
        "earthLangFunction": [
            comments("earthLangFunction"),
            {
                token: "constant.language",
                regex: "true|false"
            },
            {
                token: "string",
                regex: '"(?=.)',
                next: "qqstringType"
            },
            {
                token: "string",
                regex: "'(?=.)",
                next: "qstringType"
            },
            {
                token: "keyword.control",
                regex: "=|\\+|\\-",
            },
            {
                token: "empty",
                regex: ",",
                next: "pop"
            },
            {
                token: "empty",
                regex: "\\)",
                next: "startjs"
            },
            {
                token: "constant.language",
                regex: "i\\b",
            }, {
                token : "support.class", // hexadecimal, octal and binary
                regex : /0(?:[xX][0-9a-fA-F]+|[oO][0-7]+|[bB][01]+)\b/
            }, {
                token : "support.class", // decimal integers and floats
                regex : /(?:\d\d*(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+\b)?/
            },
            {
                token: "variable.parameter",
                regex: "\\w+\\b",
            }
        ],
        "earthLangConst" : [
            comments("earthLangConst"),
            {
                token: "entity.function.name",
                regex: "\\w+\\b"
            },
            {
                token: "empty",
                regex: "=",
                next: "startjs"
            }
        ],
        "earthLangCustomTransformParams" : [
            comments("earthLangCustomTransformParams"),
            {
                token: "storage.type",
                regex: "(?:number|complex|bool|string|array|object|function)"
            },
            {
                token: "entity.function.name",
                regex: "\\w+\\b"
            },
            {
                token: "empty",
                regex: "\\(",
            },
            {
                token: "empty",
                regex: "\\)",
            },
            {
                token: "empty",
                regex: ":",
                next: "earthLangTransform"
            },
            {
                token : "support.class", // hexadecimal, octal and binary
                regex : /0(?:[xX][0-9a-fA-F]+|[oO][0-7]+|[bB][01]+)\b/
            },
            {
                token : "support.class", // decimal integers and floats
                regex : /(?:\d\d*(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+\b)?/
            }
        ],
        "no_regex" : [
            DocCommentHighlightRules.getStartRule("doc-start"),
            comments("no_regex"),
            {
                token: "support.class.italic",
                regex: "body|camera|shader",
                push: "start",
                next: "earthLangTransform"
            },
            {
                token: "keyword.control",
                regex: "transform\\b",
                push: "start",
                next: "earthLangCustomTransformParams"
            },
            {
                token: "keyword.control",
                regex: "const\\b",
                push: "start",
                next: "earthLangConst"
            },
            {
                token: "storage.type",
                regex: "number|complex|bool|string|array|object|function",
                next: "earthLangFunction"
            },
            {
                token: "keyword.control",
                regex: "\-\>",
                next: "no_regex"
            },
            {
                token : "string",
                regex : "'(?=.)",
                next  : "qstring"
            },
            {
                token : "string",
                regex : '"(?=.)',
                next  : "qqstring"
            },
            {
                token : "constant.numeric", // hexadecimal, octal and binary
                regex : /0(?:[xX][0-9a-fA-F]+|[oO][0-7]+|[bB][01]+)\b/
            }, {
                token : "constant.numeric", // decimal integers and floats
                regex : /(?:\d\d*(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+\b)?/
            }, {
                // Sound.prototype.play =
                token : [
                    "storage.type", "punctuation.operator", "support.function",
                    "punctuation.operator", "entity.name.function", "text","keyword.operator"
                ],
                regex : "(" + identifierRe + ")(\\.)(prototype)(\\.)(" + identifierRe +")(\\s*)(=)",
                next: "function_arguments"
            }, {
                // Sound.play = function() {  }
                token : [
                    "storage.type", "punctuation.operator", "entity.name.function", "text",
                    "keyword.operator", "text", "storage.type", "text", "paren.lparen"
                ],
                regex : "(" + identifierRe + ")(\\.)(" + identifierRe +")(\\s*)(=)(\\s*)(function)(\\s*)(\\()",
                next: "function_arguments"
            }, {
                // play = function() {  }
                token : [
                    "entity.name.function", "text", "keyword.operator", "text", "storage.type",
                    "text", "paren.lparen"
                ],
                regex : "(" + identifierRe +")(\\s*)(=)(\\s*)(function)(\\s*)(\\()",
                next: "function_arguments"
            }, {
                // Sound.play = function play() {  }
                token : [
                    "storage.type", "punctuation.operator", "entity.name.function", "text",
                    "keyword.operator", "text",
                    "storage.type", "text", "entity.name.function", "text", "paren.lparen"
                ],
                regex : "(" + identifierRe + ")(\\.)(" + identifierRe +")(\\s*)(=)(\\s*)(function)(\\s+)(\\w+)(\\s*)(\\()",
                next: "function_arguments"
            }, {
                // function myFunc(arg) { }
                token : [
                    "storage.type", "text", "entity.name.function", "text", "paren.lparen"
                ],
                regex : "(function)(\\s+)(" + identifierRe + ")(\\s*)(\\()",
                next: "function_arguments"
            }, {
                // foobar: function() { }
                token : [
                    "entity.name.function", "text", "punctuation.operator",
                    "text", "storage.type", "text", "paren.lparen"
                ],
                regex : "(" + identifierRe + ")(\\s*)(:)(\\s*)(function)(\\s*)(\\()",
                next: "function_arguments"
            }, {
                // : function() { } (this is for issues with 'foo': function() { })
                token : [
                    "text", "text", "storage.type", "text", "paren.lparen"
                ],
                regex : "(:)(\\s*)(function)(\\s*)(\\()",
                next: "function_arguments"
            }, {
                // from "module-path" (this is the only case where 'from' should be a keyword)
                token : "keyword",
                regex : "from(?=\\s*('|\"))"
            }, {
                token : "keyword",
                regex : "(?:" + kwBeforeRe + ")\\b",
                next : "startjs"
            }, {
                token : ["support.constant"],
                regex : /that\b/
            }, {
                token : ["storage.type", "punctuation.operator", "support.function.firebug"],
                regex : /(console)(\.)(warn|info|log|error|time|trace|timeEnd|assert)\b/
            }, {
                token: "entity.name.function",
                regex: "^\\w+\\b",
                next: "earthLangFunctionType"
            }, {
                token : keywordMapper,
                regex : identifierRe
            }, {
                token : "punctuation.operator",
                regex : /[.](?![.])/,
                next  : "property"
            }, {
                token : "storage.type",
                regex : /=>/,
                next  : "startjs"
            }, {
                token : "keyword.control",
                regex : /--|\+\+|\.{3}|===|==|=|!=|!==|<+=?|>+=?|!|&&|\|\||\?:|[!$%&*+\-~\/^]=?/,
                next  : "startjs"
            }, {
                token : "punctuation.operator",
                regex : /[?:,;.]/,
                next  : "startjs"
            }, {
                token : "paren.lparen",
                regex : /[\[({]/,
                next  : "startjs"
            }, {
                token : "paren.rparen",
                regex : /[\])}]/
            }, {
                token: "comment",
                regex: /^#!.*$/
            },
        ],
        property: [{
                token : "text",
                regex : "\\s+"
            }, {
                // Sound.play = function play() {  }
                token : [
                    "storage.type", "punctuation.operator", "entity.name.function", "text",
                    "keyword.operator", "text",
                    "storage.type", "text", "entity.name.function", "text", "paren.lparen"
                ],
                regex : "(" + identifierRe + ")(\\.)(" + identifierRe +")(\\s*)(=)(\\s*)(function)(?:(\\s+)(\\w+))?(\\s*)(\\()",
                next: "function_arguments"
            }, {
                token : "punctuation.operator",
                regex : /[.](?![.])/
            }, {
                token : "support.function",
                regex : /(s(?:h(?:ift|ow(?:Mod(?:elessDialog|alDialog)|Help))|croll(?:X|By(?:Pages|Lines)?|Y|To)?|t(?:op|rike)|i(?:n|zeToContent|debar|gnText)|ort|u(?:p|b(?:str(?:ing)?)?)|pli(?:ce|t)|e(?:nd|t(?:Re(?:sizable|questHeader)|M(?:i(?:nutes|lliseconds)|onth)|Seconds|Ho(?:tKeys|urs)|Year|Cursor|Time(?:out)?|Interval|ZOptions|Date|UTC(?:M(?:i(?:nutes|lliseconds)|onth)|Seconds|Hours|Date|FullYear)|FullYear|Active)|arch)|qrt|lice|avePreferences|mall)|h(?:ome|andleEvent)|navigate|c(?:har(?:CodeAt|At)|o(?:s|n(?:cat|textual|firm)|mpile)|eil|lear(?:Timeout|Interval)?|a(?:ptureEvents|ll)|reate(?:StyleSheet|Popup|EventObject))|t(?:o(?:GMTString|S(?:tring|ource)|U(?:TCString|pperCase)|Lo(?:caleString|werCase))|est|a(?:n|int(?:Enabled)?))|i(?:s(?:NaN|Finite)|ndexOf|talics)|d(?:isableExternalCapture|ump|etachEvent)|u(?:n(?:shift|taint|escape|watch)|pdateCommands)|j(?:oin|avaEnabled)|p(?:o(?:p|w)|ush|lugins.refresh|a(?:ddings|rse(?:Int|Float)?)|r(?:int|ompt|eference))|e(?:scape|nableExternalCapture|val|lementFromPoint|x(?:p|ec(?:Script|Command)?))|valueOf|UTC|queryCommand(?:State|Indeterm|Enabled|Value)|f(?:i(?:nd|le(?:ModifiedDate|Size|CreatedDate|UpdatedDate)|xed)|o(?:nt(?:size|color)|rward)|loor|romCharCode)|watch|l(?:ink|o(?:ad|g)|astIndexOf)|a(?:sin|nchor|cos|t(?:tachEvent|ob|an(?:2)?)|pply|lert|b(?:s|ort))|r(?:ou(?:nd|teEvents)|e(?:size(?:By|To)|calc|turnValue|place|verse|l(?:oad|ease(?:Capture|Events)))|andom)|g(?:o|et(?:ResponseHeader|M(?:i(?:nutes|lliseconds)|onth)|Se(?:conds|lection)|Hours|Year|Time(?:zoneOffset)?|Da(?:y|te)|UTC(?:M(?:i(?:nutes|lliseconds)|onth)|Seconds|Hours|Da(?:y|te)|FullYear)|FullYear|A(?:ttention|llResponseHeaders)))|m(?:in|ove(?:B(?:y|elow)|To(?:Absolute)?|Above)|ergeAttributes|a(?:tch|rgins|x))|b(?:toa|ig|o(?:ld|rderWidths)|link|ack))\b(?=\()/
            }, {
                token : "support.function.dom",
                regex : /(s(?:ub(?:stringData|mit)|plitText|e(?:t(?:NamedItem|Attribute(?:Node)?)|lect))|has(?:ChildNodes|Feature)|namedItem|c(?:l(?:ick|o(?:se|neNode))|reate(?:C(?:omment|DATASection|aption)|T(?:Head|extNode|Foot)|DocumentFragment|ProcessingInstruction|E(?:ntityReference|lement)|Attribute))|tabIndex|i(?:nsert(?:Row|Before|Cell|Data)|tem)|open|delete(?:Row|C(?:ell|aption)|T(?:Head|Foot)|Data)|focus|write(?:ln)?|a(?:dd|ppend(?:Child|Data))|re(?:set|place(?:Child|Data)|move(?:NamedItem|Child|Attribute(?:Node)?)?)|get(?:NamedItem|Element(?:sBy(?:Name|TagName|ClassName)|ById)|Attribute(?:Node)?)|blur)\b(?=\()/
            }, {
                token :  "support.constant",
                regex : /(s(?:ystemLanguage|cr(?:ipts|ollbars|een(?:X|Y|Top|Left))|t(?:yle(?:Sheets)?|atus(?:Text|bar)?)|ibling(?:Below|Above)|ource|uffixes|e(?:curity(?:Policy)?|l(?:ection|f)))|h(?:istory|ost(?:name)?|as(?:h|Focus))|y|X(?:MLDocument|SLDocument)|n(?:ext|ame(?:space(?:s|URI)|Prop))|M(?:IN_VALUE|AX_VALUE)|c(?:haracterSet|o(?:n(?:structor|trollers)|okieEnabled|lorDepth|mp(?:onents|lete))|urrent|puClass|l(?:i(?:p(?:boardData)?|entInformation)|osed|asses)|alle(?:e|r)|rypto)|t(?:o(?:olbar|p)|ext(?:Transform|Indent|Decoration|Align)|ags)|SQRT(?:1_2|2)|i(?:n(?:ner(?:Height|Width)|put)|ds|gnoreCase)|zIndex|o(?:scpu|n(?:readystatechange|Line)|uter(?:Height|Width)|p(?:sProfile|ener)|ffscreenBuffering)|NEGATIVE_INFINITY|d(?:i(?:splay|alog(?:Height|Top|Width|Left|Arguments)|rectories)|e(?:scription|fault(?:Status|Ch(?:ecked|arset)|View)))|u(?:ser(?:Profile|Language|Agent)|n(?:iqueID|defined)|pdateInterval)|_content|p(?:ixelDepth|ort|ersonalbar|kcs11|l(?:ugins|atform)|a(?:thname|dding(?:Right|Bottom|Top|Left)|rent(?:Window|Layer)?|ge(?:X(?:Offset)?|Y(?:Offset)?))|r(?:o(?:to(?:col|type)|duct(?:Sub)?|mpter)|e(?:vious|fix)))|e(?:n(?:coding|abledPlugin)|x(?:ternal|pando)|mbeds)|v(?:isibility|endor(?:Sub)?|Linkcolor)|URLUnencoded|P(?:I|OSITIVE_INFINITY)|f(?:ilename|o(?:nt(?:Size|Family|Weight)|rmName)|rame(?:s|Element)|gColor)|E|whiteSpace|l(?:i(?:stStyleType|n(?:eHeight|kColor))|o(?:ca(?:tion(?:bar)?|lName)|wsrc)|e(?:ngth|ft(?:Context)?)|a(?:st(?:M(?:odified|atch)|Index|Paren)|yer(?:s|X)|nguage))|a(?:pp(?:MinorVersion|Name|Co(?:deName|re)|Version)|vail(?:Height|Top|Width|Left)|ll|r(?:ity|guments)|Linkcolor|bove)|r(?:ight(?:Context)?|e(?:sponse(?:XML|Text)|adyState))|global|x|m(?:imeTypes|ultiline|enubar|argin(?:Right|Bottom|Top|Left))|L(?:N(?:10|2)|OG(?:10E|2E))|b(?:o(?:ttom|rder(?:Width|RightWidth|BottomWidth|Style|Color|TopWidth|LeftWidth))|ufferDepth|elow|ackground(?:Color|Image)))\b/
            }, {
                token : "identifier",
                regex : identifierRe
            }, {
                regex: "",
                token: "empty",
                next: "no_regex"
            }
        ],
        // regular expressions are only allowed after certain tokens. This
        // makes sure we don't mix up regexps with the divison operator
        "pushedColon": [
            {
                token: "empty",
                regex: ":",
                next: "earthLangTransform"
            }
        ],

        "start": [
            comments("start"),
            {
                token: "support.class.italic",
                regex: "body|camera|shader",
                push: "start",
                next: "pushedColon"
            },
            {
                token: "storage.type",
                regex: "number|complex|bool|string|array|object|function",
                next: "earthLangFunction"
            },
            {
                token: "keyword.control",
                regex: "transform\\b",
                push: "start",
                next: "earthLangCustomTransformParams"
            },
            {
                token: "keyword.control",
                regex: "const\\b",
                push: "start",
                next: "earthLangConst"
            },
            {
                token: "entity.name.function",
                regex: "\\w+\\b",
                next: "earthLangFunctionType"
            }
        ],
        "startjs":[
            DocCommentHighlightRules.getStartRule("doc-start"),
            comments("startjs"),
            {
                token: "string.regexp",
                regex: "\\/",
                next: "regex"
            }, {
                token : "text",
                regex : "\\s+|^$",
                next : "startjs"
            }, {
                // immediately return to the start mode without matching
                // anything
                token: "empty",
                regex: "",
                next: "no_regex"
            }
        ],
        "regex": [
            {
                // escapes
                token: "regexp.keyword.operator",
                regex: "\\\\(?:u[\\da-fA-F]{4}|x[\\da-fA-F]{2}|.)"
            }, {
                // flag
                token: "string.regexp",
                regex: "/[sxngimy]*",
                next: "no_regex"
            }, {
                // invalid operators
                token : "invalid",
                regex: /\{\d+\b,?\d*\}[+*]|[+*$^?][+*]|[$^][?]|\?{3,}/
            }, {
                // operators
                token : "constant.language.escape",
                regex: /\(\?[:=!]|\)|\{\d+\b,?\d*\}|[+*]\?|[()$^+*?.]/
            }, {
                token : "constant.language.delimiter",
                regex: /\|/
            }, {
                token: "constant.language.escape",
                regex: /\[\^?/,
                next: "regex_character_class"
            }, {
                token: "empty",
                regex: "$",
                next: "no_regex"
            }, {
                defaultToken: "string.regexp"
            }
        ],
        "regex_character_class": [
            {
                token: "regexp.charclass.keyword.operator",
                regex: "\\\\(?:u[\\da-fA-F]{4}|x[\\da-fA-F]{2}|.)"
            }, {
                token: "constant.language.escape",
                regex: "]",
                next: "regex"
            }, {
                token: "constant.language.escape",
                regex: "-"
            }, {
                token: "empty",
                regex: "$",
                next: "no_regex"
            }, {
                defaultToken: "string.regexp.charachterclass"
            }
        ],
        "function_arguments": [
            {
                token: "variable.parameter",
                regex: identifierRe
            }, {
                token: "punctuation.operator",
                regex: "[, ]+"
            }, {
                token: "punctuation.operator",
                regex: "$"
            }, {
                token: "empty",
                regex: "",
                next: "no_regex"
            }
        ],
        "qqstring" : [
            {
                token : "constant.language.escape",
                regex : escapedRe
            }, {
                token : "string",
                regex : "\\\\$",
                consumeLineEnd  : true
            }, {
                token : "string",
                regex : '"|$',
                next  : "no_regex"
            }, {
                defaultToken: "string"
            }
        ],
        "qstring" : [
            {
                token : "constant.language.escape",
                regex : escapedRe
            }, {
                token : "string",
                regex : "\\\\$",
                consumeLineEnd  : true
            }, {
                token : "string",
                regex : "'|$",
                next  : "no_regex"
            }, {
                defaultToken: "string"
            }
        ]
    };


    if (!options || !options.noES6) {
        this.$rules.no_regex.unshift({
            regex: "[{}()]", onMatch: function(val, state, stack) {
                this.next = (val == "{" || val == "(") ? this.nextState : "";
                if ((val == "{" || val == "(") && stack.length) {
                    stack.unshift("js", state);
                }
                else if ((val == "}" || val == ")") && stack.length) {
                    stack.shift();
                    this.next = stack.shift();
                    if (this.next.indexOf("string") != -1 || this.next.indexOf("jsx") != -1)
                        return "paren.quasi.end";
                }
                return (val == "{" || val == "(") ? "paren.lparen" : "paren.rparen";
            },
            nextState: "startjs"
        }, {
            token : "string.quasi.start",
            regex : /`/,
            push  : [{
                token : "constant.language.escape",
                regex : escapedRe
            }, {
                token : "paren.quasi.start",
                regex : /\${/,
                push  : "startjs"
            }, {
                token : "string.quasi.end",
                regex : /`/,
                next  : "pop"
            }, {
                defaultToken: "string.quasi"
            }]
        });

        if (!options || options.jsx != false)
            JSX.call(this);
    }

    this.embedRules(DocCommentHighlightRules, "doc-",
        [ DocCommentHighlightRules.getEndRule("no_regex") ]);

    this.normalizeRules();
};

oop.inherits(EarthLangHighlightRules, TextHighlightRules);

function JSX() {
    var tagRegex = identifierRe.replace("\\d", "\\d\\-");
    var jsxTag = {
        onMatch : function(val, state, stack) {
            var offset = val.charAt(1) == "/" ? 2 : 1;
            if (offset == 1) {
                if (state != this.nextState)
                    stack.unshift(this.next, this.nextState, 0);
                else
                    stack.unshift(this.next);
                stack[2]++;
            } else if (offset == 2) {
                if (state == this.nextState) {
                    stack[1]--;
                    if (!stack[1] || stack[1] < 0) {
                        stack.shift();
                        stack.shift();
                    }
                }
            }
            return [{
                type: "meta.tag.punctuation." + (offset == 1 ? "" : "end-") + "tag-open.xml",
                value: val.slice(0, offset)
            }, {
                type: "meta.tag.tag-name.xml",
                value: val.substr(offset)
            }];
        },
        regex : "</?" + tagRegex + "",
        next: "jsxAttributes",
        nextState: "jsx"
    };
    this.$rules.start.unshift(jsxTag);
    var jsxJsRule = {
        regex: "{",
        token: "paren.quasi.start",
        push: "startjs"
    };
    this.$rules.jsx = [
        jsxJsRule,
        jsxTag,
        {include : "reference"},
        {defaultToken: "string"}
    ];
    this.$rules.jsxAttributes = [{
        token : "meta.tag.punctuation.tag-close.xml",
        regex : "/?>",
        onMatch : function(value, currentState, stack) {
            if (currentState == stack[0])
                stack.shift();
            if (value.length == 2) {
                if (stack[0] == this.nextState)
                    stack[1]--;
                if (!stack[1] || stack[1] < 0) {
                    stack.splice(0, 2);
                }
            }
            this.next = stack[0] || "startjs";
            return [{type: this.token, value: value}];
        },
        nextState: "jsx"
    },
    jsxJsRule,
    comments("jsxAttributes"),
    {
        token : "entity.other.attribute-name.xml",
        regex : tagRegex
    }, {
        token : "keyword.operator.attribute-equals.xml",
        regex : "="
    }, {
        token : "text.tag-whitespace.xml",
        regex : "\\s+"
    }, {
        token : "string.attribute-value.xml",
        regex : "'",
        stateName : "jsx_attr_q",
        push : [
            {token : "string.attribute-value.xml", regex: "'", next: "pop"},
            {include : "reference"},
            {defaultToken : "string.attribute-value.xml"}
        ]
    }, {
        token : "string.attribute-value.xml",
        regex : '"',
        stateName : "jsx_attr_qq",
        push : [
            {token : "string.attribute-value.xml", regex: '"', next: "pop"},
            {include : "reference"},
            {defaultToken : "string.attribute-value.xml"}
        ]
    },
    jsxTag
    ];
    this.$rules.reference = [{
        token : "constant.language.escape.reference.xml",
        regex : "(?:&#[0-9]+;)|(?:&#x[0-9a-fA-F]+;)|(?:&[a-zA-Z0-9_:\\.-]+;)"
    }];
}

function comments(next) {
    return [
        {
            token : "comment", // multi line comment
            regex : /\/\*/,
            next: [
                DocCommentHighlightRules.getTagRule(),
                {token : "comment", regex : "\\*\\/", next : next || "pop"},
                {defaultToken : "comment", caseInsensitive: true}
            ]
        }, {
            token : "comment",
            regex : "\\/\\/",
            next: [
                DocCommentHighlightRules.getTagRule(),
                {token : "comment", regex : "$|^", next : next || "pop"},
                {defaultToken : "comment", caseInsensitive: true}
            ]
        }
    ];
}
/* ------------------------ END ------------------------------ */
oop.inherits(EarthLangHighlightRules, TextHighlightRules);
exports.EarthLangHighlightRules = EarthLangHighlightRules;

var oop = require("../lib/oop");
var TextMode = require("./text").Mode;
var CstyleBehaviour = require("./behaviour/cstyle").CstyleBehaviour;
var CStyleFoldMode = require("./folding/cstyle").FoldMode;

var Mode = function() {
    this.HighlightRules = EarthLangHighlightRules;
    this.$behaviour = new CstyleBehaviour();
    this.foldingRules = new CStyleFoldMode();
};
oop.inherits(Mode, TextMode);

(function() {
    this.$id = "ace/mode/3arthLang";
}).call(Mode.prototype);

exports.Mode = Mode
});
