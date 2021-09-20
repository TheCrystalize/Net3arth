function parseVal(txt) {
  let v = 0;
  try {
    v = eval(txt);
    if(typeof v !== 'number') {
      throw 'NaN';
    }
  } catch (e) {
    v = parseFloat(txt);
  }
  return isNaN(v) ? 0 : v;
}

/*to*/
function complexToString(z) {
  if(z.re === 0) {
    return z.im + 'i';
  }
  if(z.im === 1) {
    return `${z.re}+i`;
  }
  if(z.im === -1) {
    return `${z.re}-i`;
  }
  if(z.im < 0) {
    return `${z.re}-${-z.im}i`;
  }
  return `${z.re}+${z.im}i`;
}

function paramsToString(params, tab) {
  let ans = '';
  for(let i = 0; i < params.length; i++) {
    switch (typeof params[i]) {
      case 'number':
      case 'string':
      case 'bool':
        ans += params[i];
        break;
      case 'object':
        if(Array.isArray(params[i])){
          let ans = '[\n'+tab;
          for(let j in params[i]) {
            ans+= '  '+paramsToString([params[i][j]], tab+'  ')+',\n'+tab;
          }
          return ans+']';
          break;
        }
        if(Object.keys(params[i]).join(',') === 're,im') {
          ans += complexToString(params[i]);
        } else if(Object.keys(params[i]).join(',') === 're,im,n') {
          if(params[i].n) {
            ans += `${params[i].re}`;
          } else {
            ans += complexToString(params[i]);
          }
        } else if(Object.keys(params[i]).join(',') === 'red,green,blue') {
          ans += `colorRGB(${params[i].red}, ${params[i].green}, ${params[i].blue})`;
        } else if(Object.keys(params[i]).join(',') === 'h,s,l') {
          ans += `colorHSL(${params[i].h}, ${params[i].s}, ${params[i].l})`;
        } else {
          ans += JSON.stringify(params[i]);
        }
        break;
      default:
        throw 'bad parameter';
    }
    if(i < params.length - 1) {
      ans += ', ';
    }
  }
  return ans;
}

function switchToString(data, tab) {
  let ans = '\n' + tab + 'choose{\n';
  for(let i = 0; i < data.length; i++) {
    ans += tab + '  ' + data[i][0] + ': ' + loopToString(data[i][1], tab + '    ') + ';\n';
  }
  return ans + tab + '}';
}

function loopToString(data, tab = '') {
  if(typeof data[0] === 'string') {
    return `${data[0]}(${paramsToString(data[1], tab)})`;
  }
  if(!data[0]) {
    return '';
  }
  if(typeof data[0][0] === 'number') {
    return switchToString(data, tab);
  }

  let ans = data.length > 1 && typeof data[1][0] === 'object' ? '\n' + tab : '';
  for(let i = 0; i < data.length; i++) {
    ans += loopToString(data[i], tab) + (i < data.length - 1 ? '\n' + tab + '-> ' : '');
  }
  return ans;
}

function customTransformsToString(transforms) {
  let ans = '';
  for(let f in transforms) {
    ans += transforms[f].name.replace(/ /g, '_') + '(' +
      transforms[f].params.map(p => p.type + ' ' + p.name + ' = ' + p.default).
    join(',') + ') {\n  ' +
      transforms[f].code.replace(/consolelog/g, 'console.log').replace(/consoleclear/g, 'console.clear') + '}\n\n';
  }
  return ans;
}

function to3arthLang(data) {
  let customFunctions = customTransformsToString(data.customFunctions);
  let body = loopToString(data.body);
  let camera = loopToString(data.camera);
  let shader = loopToString(data.shader);
  let code = customFunctions;
  if(body.length > 0) {
    code += `body: ${body};\n\n`;
  }
  if(camera.length > 0) {
    code += `camera: ${camera};\n\n`;
  }
  if(shader.length > 0) {
    code += `shader: ${shader};\n\n`;
  }
  return code;
}

/*from*/
function getTypeOfWord(word) {
  switch (true) {
    case (word[0] === '/' && word[2] === '/'):
      return 'comment';
    case (word === '->'):
    case (word === '"'):
    case (word === "'"):
    case (word === '{'):
    case (word === '}'):
    case (word === '('):
    case (word === ')'):
    case (word === '['):
    case (word === ']'):
    case (word === ':'):
    case (word === ';'):
    case (word === ','):
    case (word === '+'):
    case (word === '-'):
    case (word === '='):
    case (word === 'body'):
    case (word === 'camera'):
    case (word === 'shader'):
    case (word === 'choose'):
    case (word === 'xaos'):
      return word;
    case (word === 'number'):
    case (word === 'bool'):
    case (word === 'string'):
    case (word === 'complex'):
    case (word === 'object'):
      return 'type';
    case (word === 'true'):
    case (word === 'false'):
      return 'bool';
    case (word === 'i'):
      return 'imaginary unit';
    case (word[word.length - 1] === 'i' && word[0].search(/\d/) >= 0):
      return 'imaginary';
    case (word[0].search(/\d/) >= 0):
      return 'number';
    default:
      return 'word';
  }
}

function lineToWords(line) {
  let words = [];
  let currentWord = '';
  let typeOfWord = '';
  let isCode = false;
  let startAt = 0;

  function pushWord(word, i) {
    if(word !== '') {
      words.push({
        word: word,
        at: startAt
      });
    }
    currentWord = '';
    startAt = i;
  }
  for(let i = 0; i < line.length; i++) {
    let char = line[i];
    switch (true) {
      case (char === '#'):
        return words;
      case (char === '{'):
      case (char === '}'):
      case (char === '('):
      case (char === ')'):
      case (char === '['):
      case (char === ']'):
      case (char === ':'):
      case (char === ';'):
      case (char === ','):
      case (char === '+'):
      case (char === '='):
        pushWord(currentWord, i);
        pushWord(char, i + 1);
        typeOfWord = '';
        break;
      case (char === '-'):
        pushWord(currentWord, i);
        currentWord = char;
        typeOfWord = 'minus';
        break;
      case (typeOfWord === 'minus' && char === '>'):
        pushWord('->', i + 1);
        typeOfWord = '';
        break;
      case (char === ' '):
      case (char === '\t'):
      case (char === '\n'):
      case (char === '\r'):
        pushWord(currentWord, i + 1);
        typeOfWord = '';
        break;
      case (typeOfWord === '' && char === '.' && i < line.length-1 && line[i + 1].search(/\d/) === 0):
        currentWord = '0.';
        typeOfWord = 'decimal';
        break;
      case ((typeOfWord === 'number' || typeOfWord === 'decimal') && char.search(/\d/) === 0):
        currentWord += char;
        break;
      case (typeOfWord === 'number' && char === '.'):
        currentWord += char;
        typeOfWord = 'decimal';
        break;
      case (typeOfWord === '' && char.search(/[a-zA-Z_]/) === 0):
        currentWord = char;
        typeOfWord = 'word';
        break;
      case (typeOfWord !== 'word' && char.search(/\d/) === 0):
        pushWord(currentWord, i);
        currentWord = char;
        typeOfWord = 'number';
        break;
      case ((typeOfWord === 'number' || typeOfWord === 'decimal') && char === 'i'):
        pushWord(currentWord + 'i', i + 1);
        typeOfWord = '';
        break;
      case (typeOfWord === 'word' && char.search(/\w/) === 0):
        currentWord += char;
        break;
      default:
        pushWord(currentWord, i);
        pushWord(char, i + 1);
        typeOfWord = '';
    }
  }
  return words;
}

const htmlConsole = document.getElementById("console");

function consoleclear() {
  htmlConsole.innerText = '';
}

function escapeHtml(unsafe) {
  return unsafe
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#039;");
}

function consolelog(msg, color = "white") {
  let txt = "<div style='color:" + color + "'>";
  if(typeof msg === 'string') {
    txt += escapeHtml(msg.replace(/ /g, '\u00a0'));
  } else if(typeof msg === 'object') {
    try {
      v = escapeHtml(JSON.stringify(msg, 2));
      txt += v;
    } catch (e) {
      v = '{';
      for(let a in msg) {
        v += `\n  ${a}: ${msg[a]},`;
      }
      txt += v + '\n}';
    }
  } else if(typeof msg.toString() === 'string') {
    txt += escapeHtml(msg.toString().replace(/ /g, '\u00a0'));
  } else {
    txt += escapeHtml(msg);
  }
  htmlConsole.innerHTML += txt.replace(/\n/g, '<br/>') + "</div>";
}

function General3arthError(word, lineNumber, line) {
  return msg => {
    consolelog(`${msg}\n` +
      `Line ${lineNumber+1}:${word.at}\n` +
      `${line}\n${new Array(word.at).fill(' ').join('')}^`.replace(/ /g, '\u00a0'), "red");

    editor.getSession().setAnnotations([{
      row: lineNumber,
      column: word.at,
      text: msg,
      type: "error"
    }]);

    throw `${msg}\n` +
      `Line ${lineNumber+1}:${word.at}\n` +
      `${line}\n${new Array(word.at).fill(' ').join('')}^`;
  }
}

function _3arthError(word, lineNumber, line) {
  return msg => {
    General3arthError(word, lineNumber, line)(`Expected a ${msg}. Instead got: ${typeof word.word === 'string' ? '"'+word.word+'"' : word.word}`, lineNumber, line);
  }
}

let verbose = false;

function parseEverything(code) {
  consoleclear();
  code = code.split('\n');
  let parseState = [{
    is: "top"
  }];
  let customFunctions = {};
  for(let i = 0; i < code.length; i++) {
    let words = lineToWords(code[i]);

    wordLoop:
      for(let j = 0; j < words.length; j++) {
        let word = words[j].word;
        let wordType = getTypeOfWord(word);
        let newError = _3arthError(words[j], i, code[i]);
        let newGeneralError = General3arthError(words[j], i, code[i]);
        if(verbose) {
          console.log(customFunctions);
          console.log(JSON.stringify(parseState));
          console.log(parseState[0]);

          console.log(wordType + ': ' + word);
        }
        if(word === '/' && words[j + 1].word === '/') {
          j = words.length;
          continue;
        }

        if(word === '/' && words[j + 1].word === '*') {
          let start = [i, words[j].at];
          while(true) {
            j++;
            while(j >= words.length) {
              j = 0;
              i++;
              if(i >= code.length) {
                General3arthError({
                  at: start[1]
                }, start[0], code[start[0]])("unclosed comment");
              }
              words = lineToWords(code[i]);
            }
            if(words[j].word === '*' && words[j + 1] && words[j + 1].word === '/') {
              j++;
              continue wordLoop;
            }
          }
        }

        function endDefaultParam() {
          switch (wordType) {
            case ',':
              var param = {
                name: parseState[1].name,
                type: parseState[1].type,
                default: parseState[0].value,
              };
              parseState.shift();
              parseState.shift();
              parseState[0].params.push(param);
              break;
            case ')':
              var param = {
                name: parseState[1].name,
                type: parseState[1].type,
                default: parseState[0].value,
              };
              parseState.shift();
              parseState.shift();
              parseState[0].params.push(param);
              parseState[1].params = parseState[0].params;
              parseState.shift();
              parseState.unshift({
                is: 'js'
              });
              break;
            default:
              newError(parseState[0].is);
          }
        }

        function endParam() {
          switch (wordType) {
            case ',':
              if(!parseState[0].hasOwnProperty('value')) {
                newError(parseState[0].is);
              }
              parseState[1].params.push(parseState[0].value);
              parseState.shift();
              parseState[0].on++;
              if(parseState[0].on >= parseState[0].paramTypes.length) {
                newError('")"');
              }
              parseState.unshift({
                is: parseState[0].paramTypes[parseState[0].on].type + ' param'
              });
              break;
            case ')':
              if(!parseState[0].hasOwnProperty('value')) {
                newError(parseState[0].is);
              }
              parseState[1].params.push(parseState[0].value);
              parseState.shift();
              parseState[0].on++;
              if(parseState[0].on < parseState[0].paramTypes.length) {
                newError('","');
              }
              parseState[1].transforms.push([parseState[0].transform, parseState[0].params]);
              parseState.shift();
              parseState.unshift({
                is: 'after transform'
              });
              break;
            default:
              newError(parseState[0].is);
          }
        }

        function endWeight() {
          let terminator = parseState[0].terminator;
          if(!terminator){
            General3arthError("GET TO DA CHOPPA!");
          }
          switch (wordType) {
            case terminator:
              if(!parseState[0].hasOwnProperty('value')) {
                newError(parseState[0].is);
              }
              parseState[2].weight = parseState[0].value;
              parseState.shift();
              parseState.shift();

              parseState.unshift({
                is:"transform",
                transforms:[]
              });
              break;
            default:
              newError(parseState[0].is);
          }
        }

        function getValue() {
          let start = [i, words[j].at];
          let bracketDepth = 0;
          let parenDepth = 0;
          let squareDepth = 0;
          let lastBracketDepth = 0;
          let lastParenDepth = 0;
          let lastSquareDepth = 0;

          words = lineToWords(code[i]);

          if(words[j].word === '{') {
            bracketDepth++;
          }
          if(words[j].word === '}') {
            bracketDepth--;
          }
          if(words[j].word === '(') {
            parenDepth++;
          }
          if(words[j].word === ')') {
            parenDepth--;
          }
          if(words[j].word === '[') {
            squareDepth++;
          }
          if(words[j].word === ']') {
            squareDepth--;
          }

          if(verbose) {
            console.log(`${words[j].word} | ${parenDepth}, ${bracketDepth}, ${squareDepth}`);
          }
          do {
            lastBracketDepth = bracketDepth;
            lastParenDepth = parenDepth;
            lastSquareDepth = squareDepth;
            j++;
            while(j >= words.length) {
              j = 0;
              i++;
              if(i >= code.length) {
                General3arthError({
                  at: start[1]
                }, start[0], code[start[0]])("unexpected end of input");
              }
              words = lineToWords(code[i]);
            }
            if(words[j].word === '{') {
              bracketDepth++;
            }
            if(words[j].word === '}') {
              bracketDepth--;
            }
            if(words[j].word === '(') {
              parenDepth++;
            }
            if(words[j].word === ')') {
              parenDepth--;
            }
            if(words[j].word === '[') {
              squareDepth++;
            }
            if(words[j].word === ']') {
              squareDepth--;
            }
            if(verbose) {
              console.log(`${words[j].word} | ${parenDepth}, ${bracketDepth}, ${squareDepth}`);
            }
          } while(lastParenDepth > 0 || lastBracketDepth > 0 || lastSquareDepth > 0 || !(words[j].word === ':' || words[j].word === ',' || words[j].word === ')'))

          let jsCode = '';

          if(start[0] === i) {
            jsCode = code[i].slice(start[1], words[j].at);
          } else {
            jsCode = code[start[0]].slice(start[1]) + '\n';
            for(let l = start[0] + 1; l < i; l++) {
              jsCode += code[l] + '\n';
            }
            jsCode += code[i].slice(0, words[j].at);
          }

          jsCode = jsCode.replace(/console\.log/g, 'consolelog');
          jsCode = jsCode.replace(/console\.clear/g, 'consoleclear');

          let match = jsCode.match(/([+-]?[0-9.]+)([+-][0-9.]+)i/);
          while(match) {
            jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]})` + jsCode.slice(match.index + match[0].length);
            match = jsCode.match(/([+-]?[0-9.]+)([+-][0-9.]+)i/);
          }


          match = jsCode.match(/([+-]?[0-9.]+)([+-])i/);
          while(match) {
            jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]}1)` + jsCode.slice(match.index + match[0].length);
            match = jsCode.match(/([+-]?[0-9.]+)([+-]+)i/);
          }

          match = jsCode.match(/([-]?[0-9.]+)i/);
          while(match) {
            jsCode = jsCode.slice(0, match.index) + `C(0,${match[1]})` + jsCode.slice(match.index + match[0].length);
            match = jsCode.match(/([-]?[0-9.]+)i/);
          }

          if(verbose) {
            console.log(jsCode);
          }

          let ans;
          try {
            let loadCustomFunctions = '';
            for(let f in customFunctions) {
              loadCustomFunctions += `function ${f}(${customFunctions[f].params.map(p=>`${p.name}`).join(',')}){${customFunctions[f].code}}`;
            }
            ans = eval('(_=>{' + loadCustomFunctions + 'return ' + jsCode + '})()');
          } catch (e) {
            if(parseState[1].is === 'custom transform params') {
              newGeneralError(`${parseState[1].transform} parameter error:\n` + e);
            } else {
              newGeneralError(`${parseState[1].is} value error:\n` + e);
            }
          }

          word = words[j].word;
          wordType = getTypeOfWord(word);
          newError = _3arthError(words[j], i, code[i]);
          newGeneralError = General3arthError(words[j], i, code[i]);

          parseState[0].value = ans;

          return start;
        }

        let desiredType = parseState[0].is.replace(' param', '').replace(' weight', '');

        switch (desiredType) {
          case 'number':
          case 'complex':
          case 'bool':
          case 'string':
          case 'array':
          case 'object':
            if(verbose) {
              console.log('get value:');
            }
            let start = getValue();
            let typeError = _3arthError({
              word: parseState[0].value,
              at: start[1]
            }, start[0], code[start[0]]);
            if(typeof parseState[0].value === 'object') {
              try {
                typeError = _3arthError({
                  word: JSON.stringify(parseState[0].value),
                  at: start[1]
                }, start[0], code[start[0]]);
                if(Object.keys(parseState[0].value).join(',') === 're,im') {
                  typeError = _3arthError({
                    word: parseState[0].value.re + (parseState[0].value.im > 0 ? '+' : '') + parseState[0].value.im + 'i',
                    at: start[1]
                  }, start[0], code[start[0]]);
                }
              } catch (e) {}
            }
            if(verbose) {
              console.log(parseState[0].value);
            }
            switch (typeof parseState[0].value) {
              case 'number':
                if(desiredType !== 'number') {
                  if(desiredType === 'complex') {
                    parseState[0].value = {
                      re: parseState[0].value,
                      im: 0,
                      n: true
                    };
                  } else {
                    typeError(desiredType);
                  }
                }
                break;
              case 'object':
                if(desiredType !== 'complex') {
                  if(desiredType !== 'array' || !Array.isArray(parseState[0].value)) {
                    if(desiredType !== 'object') {
                      typeError(desiredType);
                    }
                  }
                } else if(Object.keys(parseState[0].value).join(',') === 're,im') {
                  if(typeof parseState[0].value.re !== 'number' || isNaN(parseState[0].value.re)) {
                    General3arthError({
                      word: parseState[0].value,
                      at: start[1]
                    }, start[0], code[start[0]])(`Complex value's real part is undefined: ${parseState[0].value.re}`);
                  }
                  if(typeof parseState[0].value.im !== 'number' || isNaN(parseState[0].value.im)) {
                    General3arthError({
                      word: parseState[0].value,
                      at: start[1]
                    }, start[0], code[start[0]])(`Complex value's imaginary part is undefined: ${parseState[0].value.im}`);
                  }
                } else {
                  typeError(desiredType);
                }
                break;
              case 'boolean':
                if(desiredType !== 'bool') {
                  typeError(desiredType);
                }
                break;
              case 'string':
                if(desiredType !== 'string') {
                  typeError(desiredType);
                }
                break;
              default:
                typeError(desiredType)
            }
            break;
        }

        switch (parseState[0].is) {
          case ':':
          case '=':
          case '(':
          case ')':
          case '{':
          case '}':
            if(wordType === parseState[0].is) {
              parseState.shift();
            } else {
              newError('"' + parseState[0].is + '"');
            }
            break;
          case 'top':
            switch (wordType) {
              case 'body':
                parseState.unshift({
                  is: 'body'
                });
                parseState.unshift({
                  is: 'transform',
                  transforms: []
                });
                parseState.unshift({
                  is: ':'
                });
                break;
              case 'camera':
                parseState.unshift({
                  is: 'camera'
                });
                parseState.unshift({
                  is: 'transform',
                  transforms: []
                });
                parseState.unshift({
                  is: ':'
                });
                break;
              case 'shader':
                parseState.unshift({
                  is: 'shader'
                });
                parseState.unshift({
                  is: 'transform',
                  transforms: []
                });
                parseState.unshift({
                  is: ':'
                });
                break;
              case 'word':
                parseState.unshift({
                  is: 'function',
                  name: word
                });
                break;
              default:
                newError('function, "body", "camera", or "shader"');
            }
            break;
          case 'function':
            if(wordType === '(') {
              parseState.unshift({
                is: 'function params',
                params: []
              });
            } else {
              newError('function declaration');
            }
            break;
          case 'function params':
            if(wordType === 'type') {
              parseState.unshift({
                is: 'function param',
                type: word
              });
              parseState.unshift({
                is: word
              });
              parseState.unshift({
                is: 'function param name'
              });
            } else if(wordType === ')') {
              parseState.shift();
              parseState[0].params = [];
              parseState.unshift({
                is: 'js'
              });
            } else {
              newError('type definition or ")"');
            }
            break;
          case 'function param name':
            if(wordType === 'word') {
              parseState.shift();
              parseState[1].name = word;
              parseState.unshift({
                is: '='
              });
            } else {
              newError("parameter name");
            }
            break;
          case 'function param default':
            if(wordType === 'word') {
              parseState.shift();
              parseState[0].name = word;
              parseState.unshift({
                is: 'function param default'
              });
              parseState.unshift({
                is: '='
              });
            } else {
              newError("parameter default");
            }
            break;
          case 'number':
          case 'complex':
          case 'bool':
          case 'string':
          case 'array':
          case 'object':
            endDefaultParam();
            break;
          case 'number param':
          case 'complex param':
          case 'bool param':
          case 'string param':
          case 'array param':
          case 'object param':
            endParam();
            break;
          case 'number weight':
          case 'complex weight':
          case 'bool weight':
          case 'string weight':
          case 'array weight':
          case 'object weight':
            endWeight();
            break;
          case 'js':
            if(word === '{') {
              j++;
              if(j >= words.length) {
                j = 0;
                i++;
                if(i >= code.length) {
                  throw err;
                }
                words = lineToWords(code[i]);
              }
              let err = 'unexpected end of input; unclosed "{"\n' +
                `Line ${i}:${words[j].at}\n` +
                `${code[i]}\n${new Array(words[j].at).fill(' ').join('')}^`;
              let start = [i, words[j].at];
              let bracketDepth = 1;
              while(bracketDepth > 0) {
                j++;
                while(j >= words.length) {
                  j = 0;
                  i++;
                  if(i >= code.length) {
                    throw err;
                  }
                  words = lineToWords(code[i]);
                }
                if(words[j].word === '{') {
                  bracketDepth++;
                }
                if(words[j].word === '}') {
                  bracketDepth--;
                }
              }
              parseState.shift();

              let jsCode = '';

              if(start[0] === i) {
                jsCode = code[i].slice(start[1], words[j].at);
              } else {
                jsCode = code[start[0]].slice(start[1]) + '\n';
                for(let l = start[0] + 1; l < i; l++) {
                  jsCode += code[l] + '\n';
                }
                jsCode += code[i].slice(0, words[j].at);
              }

              jsCode = jsCode.replace(/console\.log/g, 'consolelog');
              jsCode = jsCode.replace(/console\.clear/g, 'consoleclear');

              let match = jsCode.match(/([+-]?[0-9.]+)([+-][0-9.]+)i/);
              while(match) {
                jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]})` + jsCode.slice(match.index + match[0].length);
                match = jsCode.match(/([+-]?[0-9.]+)([+-][0-9.]+)i/);
              }

              match = jsCode.match(/([+-]?[0-9.]+)([+-])i/);
              while(match) {
                jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]}1)` + jsCode.slice(match.index + match[0].length);
                match = jsCode.match(/([+-]?[0-9.]+)([+-]+)i/);
              }

              match = jsCode.match(/([-]?[0-9.]+)i/);
              while(match) {
                jsCode = jsCode.slice(0, match.index) + `C(0,${match[1]})` + jsCode.slice(match.index + match[0].length);
                match = jsCode.match(/([-]?[0-9.]+)i/);
              }

              if(verbose) {
                console.log(jsCode);
              }

              try {
                let testFunction = new Function(...parseState[0].params.map(a => a.name), jsCode);
                let ans = testFunction(...parseState[0].params.map(a => a.default));
                if(typeof ans === 'function') {
                  ans = ans({
                    re: 1,
                    im: 2
                  });
                } else {
                  // determine type of result
                }
              } catch (e) {
                consolelog(`In ${parseState[0].name}()\n` + e, "red");

                editor.getSession().setAnnotations([{
                  row: start[0] - 1,
                  column: start[1],
                  text: e,
                  type: "error"
                }]);

                throw (e);
              }

              customFunctions[parseState[0].name] = {
                name: parseState[0].name,
                params: parseState[0].params,
                code: jsCode,
              };
              parseState.shift();
            } else {
              newError('JavaScript function');
            }
            break;
          case 'transform':
            switch (wordType) {
              case 'choose':
                parseState.unshift({
                  is: 'choose items',
                  items: []
                });
                parseState.unshift({
                  is: 'weight'
                });
                parseState.unshift({
                  is: 'number weight',
                  terminator: ':'
                });
                parseState.unshift({
                  is: '{'
                });
                break;
              case 'xaos':
                parseState.unshift({
                  is: 'xaos items',
                  items: []
                });
                parseState.unshift({
                  is: '{'
                });
                break;
              case 'word':
                for(let transform in BUILT_IN_TRANSFORMS_PARAMS) {
                  if(word === transform) {
                    //console.log(`found ${transform}`);
                    if(BUILT_IN_TRANSFORMS_PARAMS[word].length === 0) {
                      parseState[0].transforms.push([word, []]);
                      parseState.unshift({
                        is: 'after transform'
                      });
                      parseState.unshift({
                        is: ')'
                      });
                      parseState.unshift({
                        is: '('
                      });
                      continue wordLoop;
                    }
                    parseState.unshift({
                      is: 'custom transform params',
                      transform: word,
                      on: 0,
                      params: [],
                      paramTypes: BUILT_IN_TRANSFORMS_PARAMS[word]
                    });
                    parseState.unshift({
                      is: parseState[0].paramTypes[0].type + ' param',
                      terminator: ')'
                    });
                    parseState.unshift({
                      is: '('
                    });
                    continue wordLoop;
                  }
                }
                for(let transform in customFunctions) {
                  if(word === transform) {
                    //console.log(`found ${transform}`);
                    if(customFunctions[word].params.length === 0) {
                      parseState[0].transforms.push([word, []]);
                      parseState.unshift({
                        is: 'after transform'
                      });
                      parseState.unshift({
                        is: ')'
                      });
                      parseState.unshift({
                        is: '('
                      });
                      continue wordLoop;
                    }
                    parseState.unshift({
                      is: 'custom transform params',
                      transform: word,
                      on: 0,
                      params: [],
                      paramTypes: customFunctions[word].params
                    });
                    parseState.unshift({
                      is: parseState[0].paramTypes[0].type + ' param'
                    });
                    parseState.unshift({
                      is: '('
                    });
                    continue wordLoop;
                  }
                }
                newGeneralError(`Function undefined: ${word}`);
                break;
              default:
                newError('transform');
            }
            break;
          case 'after transform':
            switch (wordType) {
              case '->':
                parseState.shift();
                break;
              case ';':
                switch (parseState[2].is) {
                  case 'choose items':
                    parseState[2].items.push([parseState[2].weight, parseState[1].transforms]);
                    delete parseState[2].weight;
                    parseState.shift();
                    parseState.shift();
                    break;
                  case 'transform':
                    parseState[2].transforms.push(parseState[1].transforms);
                    parseState.shift();
                    parseState.shift();
                    break;
                  case 'body':
                    parseState[3].body = parseState[1].transforms;
                    parseState.shift();
                    parseState.shift();
                    parseState.shift();
                    break;
                  case 'camera':
                    parseState[3].camera = parseState[1].transforms;
                    parseState.shift();
                    parseState.shift();
                    parseState.shift();
                    break;
                  case 'shader':
                    parseState[3].shader = parseState[1].transforms;
                    parseState.shift();
                    parseState.shift();
                    parseState.shift();
                    break;
                  default:
                    newGeneralError(`Unhandeled state transition from "${parseState[1].is}" to "${parseState[2].is}"`, parseState);
                }
                break;
              default:
                newError('";" or "->"');
            }
            break;
          case 'choose items':
            switch (wordType) {
              case 'number':
                parseState[0].weight = parseFloat(word);
                parseState.unshift({
                  is: 'transform',
                  transforms: []
                });
                parseState.unshift({
                  is: ':'
                });
                break;
              case '}':
                switch (parseState[1].is) {
                  case 'choose items':
                    parseState[1].items.push([parseState[1].weight, parseState[0].items]);
                    delete parseState[1].weight;
                    parseState.shift();
                    break;
                  case 'transform':
                    parseState[1].transforms.push(parseState[0].items);
                    parseState.shift();
                    parseState.unshift({
                      is: 'after transform'
                    });
                    break;
                  default:
                    newGeneralError(`Unhandeled state transition from "${parseState[0].is}" to "${parseState[1].is}"`, parseState);
                }
                break;
              default:
                newError('weight value or "}"');
            }
            break;
          default:
            newGeneralError(`Unhandeled state: ${parseState[0].is}`, parseState);
        }
      }
  }

  if(parseState.length === 1) {
    delete parseState[0].is;

    if(!parseState[0].hasOwnProperty('body')) {
      parseState[0].body = [];
    }

    if(!parseState[0].hasOwnProperty('camera')) {
      parseState[0].camera = [];
    }

    if(!parseState[0].hasOwnProperty('shader')) {
      parseState[0].shader = [];
    }

    parseState[0].customFunctions = customFunctions;

    if(verbose) {
      console.log('-- compiled --');
    }
    return parseState[0];
  }
  General3arthError({
    at: code[code.length - 2].length - 1
  }, code.length - 1, code[code.length - 2])('unexpected end of input');
}

let autoFormat = false;

function autoFormatCode(code) {
  if(autoFormat) {
    let pos = editor.getCursorPosition();

    editor.setValue(to3arthLang(stuffToDo));
    editor.moveCursorToPosition(pos);
    editor.clearSelection();
  } else {
    editor.getSession().clearAnnotations();
  }
}

function compile3arthLang(code) {
  try {
    stuffToDo = parseEverything(code);
    autoFormatCode();
    consolelog("Finished compiling!", "limegreen");
  } catch (e) {
    console.error(e);
  }
}

let compileButton = document.getElementById('compile');

function run3arthLang(code) {
  compileButton.innerText = 'Compile';
  try {
    stuffToDo = parseEverything(code);
    autoFormatCode();
    consolelog("Finished compiling!", "limegreen");

    compileButton.innerText = 'Stop';

    refreshRender();

  } catch (e) {
    console.error(e);
  }
}

let stuffToDo;

getFile("src/default.3arth", r => {
  let pos = editor.getCursorPosition();
  editor.setValue(r);
  editor.moveCursorToPosition(pos);
  editor.clearSelection();
  //run3arthLang(r);
});
