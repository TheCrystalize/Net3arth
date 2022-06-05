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
        if(Array.isArray(params[i])) {
          let ans = '[\n' + tab;
          for(let j in params[i]) {
            ans += '  ' + paramsToString([params[i][j]], tab + '  ') + ',\n' + tab;
          }
          return ans + ']';
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

function chooseToString(data, tab) {
  let ans = 'choose{\n';
  for(let i = 0; i < data.length; i++) {
    ans += tab + '  ' + data[i][0] + ': ' + loopToString(data[i][1], tab + '    ') + ';\n';
  }
  return ans + tab + '}';
}

function xaosToString(data, tab) {
  let ans = 'xaos{\n';
  let lastEO;
  let lastAR;
  for(let i = 0; i < data.length; i++) {
    let eo = (data[i][1][0] ? (data[i][1][1] ? 'eo' : 'e') : (data[i][1][1] ? 'o' : '_'));
    let ar = '[' + paramsToString(data[i][2]) + ']';
    ans += tab + '  1:' + (eo === lastEO ? '' : eo) + ':' +
      (ar === lastAR ? '' : ar) + ':\n' + tab + '    ' +
      loopToString(data[i][3], tab + '    ') + ';\n';
    lastEO = eo;
    lastAR = ar;
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
    switch (data[0].length) {
      case 2:
        return chooseToString(data, tab);
      case 4:
        return xaosToString(data, tab);
    }
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
      transforms[f].params.map(p => p.type + ' ' + p.name).
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
    code += `body:\n${body};\n\n`;
  }
  if(camera.length > 0) {
    code += `camera:\n${camera};\n\n`;
  }
  if(shader.length > 0) {
    code += `shader:\n${shader};\n\n`;
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
    case (word === 'sum'):
    case (word === 'sumColor'):
    case (word === 'product'):
    case (word === 'productColor'):
      return word;
    case (word === 'const'):
      return 'const';
    case (word === 'transform'):
      return 'transform';
    case (word === 'number'):
    case (word === 'bool'):
    case (word === 'string'):
    case (word === 'complex'):
    case (word === 'array'):
    case (word === 'object'):
    case (word === 'function'):
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
      case (typeOfWord === '' && char === '.' && i < line.length - 1 && line[i + 1].search(/\d/) === 0):
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

function validateNumberArray(array, length, error) {
  if(!Array.isArray(array)) {
    error('weight array');
  }
  if(array.length !== length) {
    error(`${length} length weight array`);
  }
  for(let i = 0; i < length; i++) {
    if(typeof array[i] !== 'number') {
      error('numbers');
    }
  }
}

function getEnterOut(word, error) {
  switch (word) {
    case '_':
      return [0, 0];
    case 'e':
      return [1, 0];
    case 'o':
      return [0, 1];
    case 'eo':
      return [1, 1];
    default:
      error('enterOut definition');
  }
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

function getPrecomputeString(preComputeStuff) {
  let preComputeString = '';
  for(let thing in preComputeStuff) {
    switch (preComputeStuff[thing].is) {
      case '_function':
        preComputeString += `function ${preComputeStuff[thing].name}(${preComputeStuff[thing].params.map(p=>`${p.name}`).join(',')}){${preComputeStuff[thing].code}}`;
        break;
      case 'const':
        preComputeString += `const ${preComputeStuff[thing].name} = ${preComputeStuff[thing].const};`;
        break;
    }
  }
  return preComputeString;
}

function parseEverything(code) {
  code = code.split('\n');
  let parseState = [{
    is: "top"
  }];
  let preComputeStuff = [];
  let consts = {};
  let customFunctions = {};
  let transformTemplates = {};
  let customTransformParams = false;
  for(let i = 0; i < code.length; i++) {
    let words = lineToWords(code[i]);

    wordLoop:
      for(let j = 0; j < words.length; j++) {
        let word = words[j].word;
        let wordType = getTypeOfWord(word);
        let newError = _3arthError(words[j], i, code[i]);
        let newGeneralError = General3arthError(words[j], i, code[i]);
        if(verbose) {
          console.log('-- loop-top --');
          console.log(consts);
          console.log(transformTemplates);
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

        function evalValue(jsCode) {
          if(verbose) {
            console.log(`~~ eval value ~~`);
            console.log(typeof jsCode);
          }
          if(typeof jsCode !== 'string') {
            if(typeof jsCode === 'object') {
              return {
                code: `(${JSON.stringify(jsCode)})`,
                isParam: false
              }
            }
            return {
              code: `(${jsCode})`,
              isParam: false
            }
          }
          jsCode = jsCode.replace(/console\.log/g, 'consolelog');
          jsCode = jsCode.replace(/console\.clear/g, 'consoleclear');

          let match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-][0-9]+[0-9.]*)i/);
          while(match) {
            jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]})` + jsCode.slice(match.index + match[0].length);
            match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-][0-9]+[0-9.]*)i/);
          }

          match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-])i/);
          while(match) {
            jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]}1)` + jsCode.slice(match.index + match[0].length);
            match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-]+)i/);
          }

          match = jsCode.match(/([-]?[0-9]+[0-9.]*)i/);
          while(match) {
            jsCode = jsCode.slice(0, match.index) + `C(0,${match[1]})` + jsCode.slice(match.index + match[0].length);
            match = jsCode.match(/([-]?[0-9]+[0-9.]*)i/);
          }

          if(verbose) {
            console.log(jsCode);
          }

          let isParam = false;
          let testCode = jsCode;
          if(customTransformParams) {
            for(let i = 0; i < customTransformParams.length; i++) {
              if(jsCode.match(new RegExp(`\\b${customTransformParams[i].name}\\b`))) {
                jsCode = jsCode.replace(new RegExp(`\\b${customTransformParams[i].name}\\b`, 'g'), `<!${customTransformParams[i].name}!>`);
                let replacement = '';
                switch (customTransformParams[i].type) {
                  case 'number':
                    replacement = '0';
                    break;
                  case 'complex':
                    replacement = '{re:0,im:0}';
                    break;
                  case 'bool':
                    replacement = 'false';
                    break;
                  case 'string':
                    replacement = '"test"';
                    break;
                  case 'array':
                    replacement = '[0]';
                    break;
                  case 'object':
                    replacement = '{re: 0, im: 0, z: 0, red: 0, green: 0, blue: 0, alpha: 0, h: 0, s: 0, l: 0}';
                    break;
                  case 'function':
                    replacement = '_=>0';
                    break;
                  default:
                    newGeneralError(`Unhandeled type test: ${customTransformParams[i].type}`);
                }
                testCode = testCode.replace(new RegExp(`\\b${customTransformParams[i].name}\\b`, 'g'), '(' + replacement + ')');
                isParam = true;
              }
            }
          }

          let preComputeString = getPrecomputeString(preComputeStuff);

          return {
            code: '(_=>{return ' + jsCode + '})()',
            testCode: '(_=>{' + preComputeString + 'return ' + testCode + '})()',
            isParam: isParam
          };
        }

        function replaceParams(params, values, transforms) {
          let ans = [];
          for(let i = 0; i < transforms.length; i++) {
            ans.push(replaceParam(params, values, transforms[i]));
          }
          return ans;
        }

        function customChoose(params, values, transform) {
          let ans = ['choose'];
          for(let i = 1; i < transform.length; i++) {
            ans.push([
              replaceParam(params, values, transform[i][0]),
              customTransform(params, values, transform[i][1])
            ]);
          }
          return ans;
        }

        function customXaos(params, values, transform) {
          let ans = ['xaos'];
          for(let i = 1; i < transform.length; i++) {
            ans.push([
              replaceParam(params, values, transform[i][0]),
              transform[i][1],
              replaceParam(params, values, transform[i][2]),
              customTransform(params, values, transform[i][3])
            ]);
          }
          return ans;
        }

        function customSum(params, values, transform) {
          let ans = ['sum'];
          for(let i = 1; i < transform.length; i++) {
            ans.push([
              replaceParam(params, values, transform[i][0]),
              transform[i][1],
              replaceParam(params, values, transform[i][2]),
              customTransform(params, values, transform[i][3])
            ]);
          }
          return ans;
        }

        function customSumColor(params, values, transform) {
          let ans = ['sumColor'];
          for(let i = 1; i < transform.length; i++) {
            ans.push([
              replaceParam(params, values, transform[i][0]),
              transform[i][1],
              replaceParam(params, values, transform[i][2]),
              customTransform(params, values, transform[i][3])
            ]);
          }
          return ans;
        }

        function customProduct(params, values, transform) {
          let ans = ['product'];
          for(let i = 1; i < transform.length; i++) {
            ans.push([
              replaceParam(params, values, transform[i][0]),
              transform[i][1],
              replaceParam(params, values, transform[i][2]),
              customTransform(params, values, transform[i][3])
            ]);
          }
          return ans;
        }

        function customProductColor(params, values, transform) {
          let ans = ['productColor'];
          for(let i = 1; i < transform.length; i++) {
            ans.push([
              replaceParam(params, values, transform[i][0]),
              transform[i][1],
              replaceParam(params, values, transform[i][2]),
              customTransform(params, values, transform[i][3])
            ]);
          }
          return ans;
        }

        function customTransform(params, values, transform) {
          if(typeof transform[0] === 'string') {
            return [transform[0], replaceParams(params, values, transform[1])];
          }
          if(!transform[0]) {
            return '';
          }
          if(typeof transform[0][0] === 'string') {
            switch (transform[0][0]) {
              case 'choose':
                return customChoose(params, values, transform);
                break;
              case 'xaos':
                return customXaos(params, values, transform);
                break;
              case 'sum':
                return customSum(params, values, transform);
                break;
              case 'sumColor':
                return customSumColor(params, values, transform);
                break;
              case 'product':
                return customProduct(params, values, transform);
                break;
              case 'productColor':
                return customProductColor(params, values, transform);
                break;
            }
          }
          if(typeof transform[0][0] === 'number') {
            console.log(transform);
            console.log(transform[0]);
            throw 'bad number';
            switch (transform[0].length) {
              case 2:
                return customChoose(params, values, transform);
              case 4:
                return customXaos(params, values, transform);
            }
          }

          let ans = [];
          for(let i = 0; i < transform.length; i++) {
            ans.push(customTransform(params, values, transform[i]));
          }
          return ans;
        }

        function replaceParam(params, values, txt) {
          if(typeof txt !== 'string') {
            return txt;
          }
          let ans = txt;
          for(let i = 0; i < values.length; i++) {
            ans = ans.replace(new RegExp('<!' + params[i].name + '!>', 'g'), values[i]);
          }

          if(ans.indexOf('<!') >= 0) {
            return ans;
          }

          return eval(evalValue(ans).testCode);
        }

        function postParams() {
          switch (parseState[0].is) {
            case '_function':
              parseState.unshift({
                is: 'js'
              });
              break;
            case 'new transform':
              customTransformParams = parseState[0].params;
              parseState.unshift({
                is: 'transform',
                transforms: []
              });
              parseState.unshift({
                is: ':'
              });
              break;
            default:
              newGeneralError(`Unhandeled state transition from "parameters" "${parseState[0].is}"`, parseState);
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
              switch (parseState[0].is) {
                case 'transform':
                  parseState.unshift({
                    is: 'after transform'
                  });
                  break;
                case 'parse custom function':
                  parseState[1].transforms.push(...customTransform(
                    transformTemplates[parseState[0].name].params,
                    parseState[0].transforms[0][1],
                    transformTemplates[parseState[0].name].transform
                  ));
                  parseState.shift();
                  parseState.unshift({
                    is: 'after transform'
                  });
                  break;
                default:
                  newGeneralError(`Unhandeled state transition: endParam -> ${parseState[0].is}`);
              }
              break;
            default:
              newError(parseState[0].is);
          }
        }

        function endWeight() {
          let terminator = parseState[0].terminator;
          if(!terminator) {
            newGeneralError("GET TO DA CHOPPA!");
          }
          switch (wordType) {
            case terminator:
              if(!parseState[0].hasOwnProperty('value')) {
                newError(parseState[0].is);
              }
              parseState[2].weight = parseState[0].value;
              switch (parseState[2].is) {
                case 'sum items':
                case 'sumColor items':
                case 'product items':
                case 'productColor items':
                case 'choose items':
                  parseState.shift();
                  parseState.shift();

                  parseState.unshift({
                    is: "transform",
                    transforms: []
                  });
                  break;
                case 'xaos items':
                  switch (parseState[0].is) {
                    case 'number weight':
                      parseState.shift();
                      parseState.shift();
                      break;
                    case 'array weight':
                      parseState.shift();
                      parseState.shift();
                      if(parseState[0].items.length === 0) {
                        validateNumberArray(parseState[0].weight, parseState[0].weight.length, newError);
                        parseState[0].size = parseState[0].weight.length;
                      } else {
                        validateNumberArray(parseState[0].weight, parseState[0].size, newError);
                      }
                      parseState.unshift({
                        is: "transform",
                        transforms: []
                      });
                      break;
                  }
                  break;
                default:
                  newGeneralError(`Unhandeled endWeight for "${parseState[2].is}" "${parseState[0].is}"`, parseState);
              }
              break;
            default:
              newError(parseState[0].is);
          }
        }

        function getValue() {
          let ans;
          let start = [i, words[j].at];
          if(word === ':') {
            switch (parseState[2].is) {
              case 'xaos items':
                if(parseState[2].items.length > 0) {
                  if(!parseState[2].hasOwnProperty('enterOut')) {
                    ans = parseState[2].items[parseState[2].items.length - 1][1];
                  } else {
                    ans = parseState[2].items[parseState[2].items.length - 1][2];
                  }
                } else {
                  newGeneralError("No ");
                }
                parseState[0].value = ans;
                return start;
                return;
              default:
                newGeneralError("Unexpected token ':'");
                return;
            }
          }
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
          } while(lastParenDepth > 0 || lastBracketDepth > 0 || lastSquareDepth > 0 || !(words[j].word === ';' || words[j].word === ':' || words[j].word === ',' || words[j].word === ')'))

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

          jsCode = evalValue(jsCode);

          try {
            if(verbose) {
              console.log(jsCode);
              console.log(`isParam: ${jsCode.isParam}`);
            }
            if(jsCode.isParam) {
              let test = eval(jsCode.testCode);
              ans = jsCode.code;
              parseState[0].testValue = test;
            } else {
              ans = eval(jsCode.testCode);
            }
          } catch (e) {
            if(verbose) {
              //console.log(JSON.stringify(parseState));
            }
            if(parseState[1].is === 'custom transform params') {
              newGeneralError(`${parseState[1].transform} parameter error:\n` + e);
            } else {
              newGeneralError(`${parseState[1].is} value error:\n` + e);
            }
          }

          if(typeof ans === 'function'){
            ans = {functionCode: ans + ''};
          }

          word = words[j].word;
          wordType = getTypeOfWord(word);
          newError = _3arthError(words[j], i, code[i]);
          newGeneralError = General3arthError(words[j], i, code[i]);

          parseState[0].value = ans;

          return start;
        }

        if(parseState[0].is.indexOf(' items') > 0 && word !== '}' && word !== ';' && !parseState[0].hasOwnProperty('weight')) {
          parseState.unshift({
            is: 'weight'
          });
          parseState.unshift({
            is: 'number weight',
            terminator: ':'
          });
        }

        let desiredType = parseState[0].is.replace(' param', '').replace(' weight', '');

        switch (desiredType) {
          case 'number':
          case 'complex':
          case 'bool':
          case 'string':
          case 'array':
          case 'object':
          case 'function':
            if(verbose) {
              console.log('get value:');
            }
            let start = getValue();

            let val = parseState[0].value;

            if(customTransformParams) {
              for(let i = 0; i < customTransformParams.length; i++) {
                if(typeof parseState[0].value !== 'string') {
                  continue;
                }
                if(parseState[0].value.indexOf('<!' + customTransformParams[i].name + '!>') >= 0) {
                  val = parseState[0].testValue;
                  i = Infinity;
                }
              }

              if(desiredType === 'complex' && typeof val === 'number') {
                if(typeof parseState[0].value === 'string' && parseState[0].value.indexOf('<!') >= 0) {
                  parseState[0].value = 'C(' + parseState[0].value + ',0)';
                } else {
                  parseState[0].value = eval('C(' + parseState[0].value + ',0)');
                }
              }
            }

            if(desiredType === 'complex' && typeof parseState[0].value === 'number') {
              parseState[0].value = {
                re: parseState[0].value,
                im: 0
              };
            }

            if(verbose) {
              console.log(customTransformParams);
              console.log(start);
              console.log(JSON.stringify(parseState));
              console.log(val);
            }
            let typeError = _3arthError({
              word: val,
              at: start[1]
            }, start[0], code[start[0]]);
            if(typeof val === 'object') {
              try {
                typeError = _3arthError({
                  word: JSON.stringify(val),
                  at: start[1]
                }, start[0], code[start[0]]);
                if(Object.keys(val).join(',') === 're,im') {
                  typeError = _3arthError({
                    word: val.re + (val.im > 0 ? '+' : '') + val.im + 'i',
                    at: start[1]
                  }, start[0], code[start[0]]);
                }
              } catch (e) {}

              if(verbose) {
                console.log(val);
              }
              switch (typeof val) {
                case 'number':
                  if(desiredType !== 'number') {
                    if(desiredType === 'complex') {
                      val = {
                        re: val,
                        im: 0,
                        n: true
                      };
                    } else {
                      typeError(desiredType);
                    }
                  }
                  break;
                case 'object':
                  if(desiredType !== 'complex' && desiredType !== 'function') {
                    if(desiredType !== 'array' || !Array.isArray(val)) {
                      if(desiredType !== 'object') {
                        typeError(desiredType);
                      }
                    }
                  }
                  else if(desiredType !== 'complex'){
                    if(!val.hasOwnProperty('functionCode')){
                      typeError(desiredType);
                    }
                  }
                  else if(Object.keys(val).join(',') === 're,im') {
                    if(typeof val.re !== 'number' || isNaN(val.re)) {
                      General3arthError({
                        word: val,
                        at: start[1]
                      }, start[0], code[start[0]])(`Complex value's real part is undefined: ${val.re}`);
                    }
                    if(typeof val.im !== 'number' || isNaN(val.im)) {
                      General3arthError({
                        word: val,
                        at: start[1]
                      }, start[0], code[start[0]])(`Complex value's imaginary part is undefined: ${val.im}`);
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
            }
            break;
        }

        if(verbose) {
          console.log('-- switch-top --');
          console.log(transformTemplates);
          console.log(customFunctions);
          console.log(JSON.stringify(parseState));
          console.log(parseState[0]);

          console.log(wordType + ': ' + word);
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
              case 'const':
                parseState.unshift({
                  is: 'new const'
                });
                parseState.unshift({
                  is: 'new name'
                });
                break;
              case 'transform':
                parseState.unshift({
                  is: 'new transform'
                });
                parseState.unshift({
                  is: 'new name'
                });
                break;
              case 'word':
                parseState.unshift({
                  is: '_function',
                  name: word
                });
                break;
              default:
                newError('function, "body", "camera", or "shader"');
            }
            break;
          case '_function':
            if(wordType === '(') {
              parseState.unshift({
                is: '_function params',
                params: []
              });
            } else {
              newError('_function declaration');
            }
            break;
          case '_function param':
            switch (wordType) {
              case ',':
                parseState.shift();
                break;
              case ')':
                parseState.shift();
                parseState[1].params = parseState[0].params;
                parseState.shift();
                postParams();
                break;
              default:
                newError('"," or ")"');
            }
            break;
          case '_function params':
            if(wordType === 'type') {
              parseState.unshift({
                is: '_function param',
                type: word
              });
              parseState.unshift({
                is: '_function param name'
              });
            } else if(wordType === ')') {
              parseState.shift();
              parseState[0].params = [];
              postParams();
            } else {
              newError('type definition or ")"');
            }
            break;
          case '_function param name':
            if(wordType === 'word') {
              parseState.shift();
              parseState[1].params.push({
                name: word,
                type: parseState[0].type
              });
            } else {
              newError("parameter name");
            }
            break;
          case 'number param':
          case 'complex param':
          case 'bool param':
          case 'string param':
          case 'array param':
          case 'object param':
          case 'function param':
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
              jsCode = jsCode.replace(/console\.log/g, 'consolelog');
              jsCode = jsCode.replace(/\bzBuffer\( *\)/g, 'zBuffer(oldBuffer, newBuffer)');
              jsCode = jsCode.replace(/\bfirstBuffer\( *\)/g, 'firstBuffer(oldBuffer, newBuffer)');
              jsCode = jsCode.replace(/\blastBuffer\( *\)/g, 'lastBuffer(oldBuffer, newBuffer)');
              jsCode = jsCode.replace(/\baverageBuffer\( *\)/g, 'averageBuffer(oldBuffer, newBuffer)');

              let match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-][0-9]+[0-9.]*)i\b/);
              while(match) {
                jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]})` + jsCode.slice(match.index + match[0].length);
                match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-][0-9]+[0-9.]*)i\b/);
              }

              match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-])i\b/);
              while(match) {
                jsCode = jsCode.slice(0, match.index) + `C(${match[1]},${match[2]}1)` + jsCode.slice(match.index + match[0].length);
                match = jsCode.match(/([+-]?[0-9]+[0-9.]*)([+-]+)i\b/);
              }

              match = jsCode.match(/([-]?[0-9]+[0-9.]*)i\b/);
              while(match) {
                jsCode = jsCode.slice(0, match.index) + `C(0,${match[1]})` + jsCode.slice(match.index + match[0].length);
                match = jsCode.match(/([-]?[0-9]+[0-9.]*)i\b/);
              }

              if(verbose) {
                console.log(jsCode);
              }

              if(parseState[0].name === 'buffer') {
                parseState[0].params = [{
                  name: "oldBuffer",
                  type: "object",
                }, {
                  name: "newBuffer",
                  type: "object"
                }];
              }

              try {
                let ans;
                {
                  let preComputeString = getPrecomputeString(preComputeStuff);
                  let testfunction = new Function(...parseState[0].params.map(a => a.name), preComputeString + jsCode);
                  let sampleParams = parseState[0].params.map(a => {
                    switch (a.type) {
                      case 'string':
                        return 'test';
                      case 'boolean':
                        return false;
                      case 'number':
                        return 0;
                      case 'array':
                        return [];
                      case 'complex':
                        return C(0, 0);
                      case 'object':
                        return {
                          red: 0, green: 0, blue: 0, h: 0, s: 0, l: 0, re: 0, im: 0, z: 0, alpha: 0
                        };
                      case 'function':
                        return _ => 0;
                      default:
                        return 0;
                    }
                  });
                  ans = testfunction(...sampleParams);
                }
                if(parseState[0].name === 'buffer') {
                  if(!ans.hasOwnProperty('red') || !ans.hasOwnProperty('green') || !ans.hasOwnProperty('blue') || !ans.hasOwnProperty('z')) {
                    throw `buffer must return at least {red, green, blue, z}\ngot:\n${typeof ans === 'object'?JSON.stringify(ans,0,2):ans}`;
                  }
                } else {
                  // TODO: determine type of result?
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
                is: '_function',
                name: parseState[0].name,
                params: parseState[0].params,
                code: jsCode,
              };
              preComputeStuff.push(customFunctions[parseState[0].name]);
              parseState.shift();
            } else {
              newError('JavaScript function');
            }
            break;
          case 'transform':
            switch (wordType) {
              case 'sum':
              case 'sumColor':
              case 'product':
              case 'productColor':
              case 'choose':
                parseState.unshift({
                  is: `${wordType} items`,
                  items: [wordType]
                });
                parseState.unshift({
                  is: '{'
                });
                break;
              case 'xaos':
                parseState.unshift({
                  is: `${wordType} items`,
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
                for(let transform in transformTemplates) {
                  if(word === transform) {
                    //console.log(`found ${transform}`);
                    if(transformTemplates[word].params.length === 0) {
                      parseState[0].transforms.push(...transformTemplates[word].transform);
                      parseState.unshift({
                        is: 'after transform',
                        transforms: []
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
                      is: 'parse custom function',
                      name: word,
                      transforms: []
                    });
                    parseState.unshift({
                      is: 'transform template params',
                      transform: word,
                      on: 0,
                      params: [],
                      paramTypes: transformTemplates[word].params
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
                newGeneralError(`function undefined: ${word}`);
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
                  case 'new transform':
                    transformTemplates[parseState[2].name] = {
                      params: parseState[2].params,
                      transform: parseState[1].transforms
                    };
                    customTransformParams = false;
                    parseState.shift();
                    parseState.shift();
                    parseState.shift();
                    break;
                  case 'sum items':
                  case 'sumColor items':
                  case 'product items':
                  case 'productColor items':
                  case 'choose items':
                    parseState[2].items.push([parseState[2].weight, parseState[1].transforms]);
                    delete parseState[2].weight;
                    parseState.shift();
                    parseState.shift();
                    break;
                  case 'xaos items':
                    parseState[2].items.push([parseState[2].mainWeight, parseState[2].enterOut, parseState[2].weight, parseState[1].transforms]);
                    delete parseState[2].mainWeight;
                    delete parseState[2].enterOut;
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
                    newGeneralError(`Unhandeled state transition after "${parseState[1].is}" to "${parseState[2].is}"`, parseState);
                }
                break;
              default:
                newError('";" or "->"');
            }
            break;
          case 'sum items':
          case 'sumColor items':
          case 'product items':
          case 'productColor items':
          case 'choose items':
            switch (wordType) {
              case '}':
                switch (parseState[1].is) {
                  case parseState[0].is:
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
          case 'xaos items':
            if(parseState[0].hasOwnProperty('weight') &&
              !parseState[0].hasOwnProperty('enterOut')) {
              parseState[0].mainWeight = parseState[0].weight;
              delete parseState[0].weight;
              let same = word === ':';
              if(same) {
                if(parseState[0].items.length === 0) {
                  newGeneralError("Unexpected token ':'");
                }
                parseState[0].enterOut = parseState[0].items[parseState[0].items.length - 1][1];
              } else {
                parseState[0].enterOut = getEnterOut(word, newError);
              }
              parseState.unshift({
                is: 'weight'
              });
              parseState.unshift({
                is: 'array weight',
                terminator: ':'
              });
              if(!same) {
                parseState.unshift({
                  is: ':'
                });
              }
              break;
            }
            switch (wordType) {
              case '}':
                if(parseState[0].size < parseState[0].items.length) {
                  newGeneralError(`Xaos weight arrays are only length ${parseState[0].size} but there ${parseState[0].items.length} states.`);
                }
                if(parseState[0].size > parseState[0].items.length) {
                  newGeneralError(`Xaos weight arrays are length ${parseState[0].size} but there are only ${parseState[0].items.length} states.`);
                }
                switch (parseState[1].is) {
                  case 'transform':
                    parseState[1].transforms.push(['xaos', ...JSON.parse(JSON.stringify(parseState[0].items)).map(row => {
                      let newWeights = row[2].map((sub, index) => parseState[0].items[index][0] * sub);
                      return [newWeights.reduce((a, b) => a + b), row[1], newWeights, row[3]];
                    })]);
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
                newError(`}`);
            }
            break;

          case 'variable param':
            let val = getValue();
            try {
              if(verbose) {
                console.log(`~~ variable ~~`);
                console.log(parseState[0].value);
                console.log(evalValue(parseState[0].value).code);
                console.log(eval(evalValue(parseState[0].value).code));
              }
              parseState[0].value = eval(evalValue(parseState[0].value).code);
            } catch (e) {
              newGeneralError('declaration error:\n' + e);
            }
            let cv = parseState[0].value;

            switch (typeof cv) {
              case 'object':
                cv = JSON.stringify(cv);
                break;
              default:
                cv = cv + '';
            }

            consts[parseState[1].name] = cv;
            preComputeStuff.push({
              is: 'const',
              name: parseState[1].name,
              const: cv
            });
            parseState.shift();
            parseState.shift();
            break;

          case 'new name':
            switch (wordType) {
              case 'word':
                parseState[1].name = word;
                parseState.shift();
                switch (parseState[0].is) {
                  case 'new transform':
                    parseState.unshift({
                      is: '_function params',
                      params: []
                    });
                    parseState.unshift({
                      is: '('
                    });
                    break;
                  case 'new const':
                    parseState.unshift({
                      is: 'variable param',
                      terminator: ';'
                    });
                    parseState.unshift({
                      is: '='
                    });
                    break;
                  default:
                    newGeneralError('Something went horribly wrong...');
                }
                break;
              default:
                newError('valid alphanumeric name');
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

    parseState[0].preCompute = preComputeStuff;

    if(verbose) {
      console.log('-- compiled --');
      console.log(JSON.stringify(parseState[0]));
    }
    return JSON.parse(JSON.stringify(parseState[0]));
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
  } catch (e) {
    console.error(e);
    runButton.innerText = 'Run';
    compileButton.innerText = 'Compile';
    throw ('Compile error');
  }
}

let compileButton = document.getElementById('compile');
let runButton = document.getElementById('run');

function resume3arthLang(code) {
  compileButton.innerText = 'Compile';
  try {
    stuffToDo = parseEverything(code);
    autoFormatCode();

    compileButton.innerText = 'Stop';

    refreshRender(false);

  } catch (e) {
    console.error(e);
    runButton.innerText = 'Run';
    compileButton.innerText = 'Compile';
    throw ('Compile error');
  }
}

function run3arthLang(code) {
  compileButton.innerText = 'Compile';
  try {
    consoleclear();
    stuffToDo = parseEverything(code);
    autoFormatCode();

    consolelog("Finished compiling!", "limegreen");

    compileButton.innerText = 'Stop';

    refreshRender();

  } catch (e) {
    console.error(e);
    runButton.innerText = 'Run';
    compileButton.innerText = 'Compile';
    throw ('Compile error');
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
