let editor = ace.edit("editor");
editor.setTheme("ace/theme/tomorrow_night");
editor.session.setMode("ace/mode/3arthLang");

editor.setOptions({
  tabSize: 2,
  fontFamily: "Fira Code",
  fontSize: "16px",
  wrap: true,
  enableBasicAutocompletion: true,
  enableLiveAutocompletion: true,
  useWorker: false
});

var customWordCompleter = {
  getCompletions: function(editor, session, pos, prefix, callback) {
    var wordList = [
      ["TEMPLATE", ["JavaScript TRANSFORM"]],
      ["BUFFER", ["JavaScript BUFFER", "averageBuffer", "firstBuffer", "lastBuffer", "zBuffer"]],
      ["CONST", ["const"]],
      ["TRANSFORM", ["transform"]],
      ["BODY", ["body"]],
      ["CAMERA", ["camera"]],
      ["SHADER", ["shader"]],
      ["XAOS", ["xaos{}"]],
      ["CHOOSE", ["choose{}"]],
      ["SUM", ["sum{}"]],
      ["SUM COLOR", ["sumColor{}"]],
      ["PRODUCT", ["product{}"]],
      ["PRODUCT COLOR", ["productColor{}"]],
      ["keyword", ["choose", "sum", "sumColor", "product", "productColor", "xaos"]],
      ["type", ["complex", "number", "string", "bool", "array", "object"]],
      ["constructor", StandardLibConstructorNamesArray],
      ["constant", StandardLibConstantsNamesArray],
      ["transform", StandardLibNamesArray.map(a=>a+'()')],
      ["helper", StandardLibHelperNamesArray],
      ["keyword", [...StandardLibNamesArray]],
    ];
    callback(null, wordList.map(function(words) {
      return words[1].map(function(word) {
        let post = word;
        switch (words[0]) {
          case "CONST":
          case "TRANSFORM":
            post += ' ';
            break;
          case "CAMERA":
          case "BODY":
          case "SHADER":
            post += ":\n";
            break;
          case "TEMPLATE":
            post = "myTransform() {\n  return z => {\n    return {\n      ...z,\n      re: z.re,\n      im: z.im\n    }\n  }\n}\n";
            break;
          case "BUFFER":
            if(word === "JavaScript BUFFER"){
              post = "buffer() {\n  return {\n    red: oldBuffer.red + newBuffer.red * newBuffer.alpha,\n    green: oldBuffer.green + newBuffer.green * newBuffer.alpha,\n    blue: oldBuffer.blue + newBuffer.blue * newBuffer.alpha,\n    z: 0,\n  }\n}\n";}
            else{
              post = `buffer() {\n  return ${word}();\n}\n`;
            }
            break;
          case "XAOS":
            post = "xaos{\n  1:eo:[1,1]:\n}";
            break;
          case "CHOOSE":
            post = "choose{\n  1:\n}";
            break;
          case "SUM":
            post = "sum{\n  1:\n}";
            break;
          case "SUM COLOR":
            post = "sumColor{\n  1:\n}";
            break;
          case "PRODUCT":
            post = "product{\n  1:\n}";
            break;
          case "PRODUCT COLOR":
            post = "productColor{\n  1:\n}";
            break;
          case "type":
            post += " ";
            break;
          case "constructor":
          case "helper":
            post += "(";
            break;
          case "transform":
            if(BUILT_IN_TRANSFORMS_PARAMS.hasOwnProperty(word.slice(0,word.length-2))) {
              wordT = word.slice(0,word.length-2);
              post = wordT;
              post += "(";
              for(let i = 0; i < BUILT_IN_TRANSFORMS_PARAMS[wordT].length; i++) {
                post += paramsToString([BUILT_IN_TRANSFORMS_PARAMS[wordT][i].default],'');
               if(i < BUILT_IN_TRANSFORMS_PARAMS[wordT].length - 1) {
                  post += ", ";
                }
              }
              post += ")";
            }
            break;
        };
        return {
          caption: word,
          value: post,
          meta: words[0],
          score: Infinity
        };
      })
    }).flat());
  }
}

editor.completers.unshift(customWordCompleter);


editor.commands.addCommand({
  name: "run",
  bindKey: {
    win: "Ctrl-R|Ctrl-Space",
    mac: "Command-R|Command-Space"
  },
  exec: function() {
    runCode();
  }
});


editor.commands.addCommand({
  name: "stop",
  bindKey: {
    win: "Ctrl-Shift-R|Ctrl-S",
    mac: "Command-Shift-R|Command-S"
  },
  exec: function() {
    stopCode();
  }
});

editor.commands.addCommand({
  name: "next tutorial",
  bindKey: {
    win: "Ctrl-/",
    mac: "Command-/"
  },
  exec: function() {
    tutorial();
  }
});

editor.commands.addCommand({
  name: "last tutorial",
  bindKey: {
    win: "Ctrl-Shift-/",
    mac: "Command-Shift-/"
  },
  exec: function() {
    tutorial(true);
  }
});
