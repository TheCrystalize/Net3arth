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
      ["TEMPLATE", ["TRANSFORM"]],
      ["BODY", ["body"]],
      ["CAMERA", ["camera"]],
      ["SHADER", ["shader"]],
      ["keyword", ["choose", "xaos"]],
      ["type", ["complex", "number", "string", "bool", "array", "object"]],
      ["constructor", StandardLibConstructorNamesArray],
      ["transform", StandardLibNamesArray.map(a=>a+'()')],
      ["helper", StandardLibHelperNamesArray],
      ["keyword", [...StandardLibNamesArray]],
    ];
    callback(null, wordList.map(function(words) {
      return words[1].map(function(word) {
        let post = word;
        switch (words[0]) {
          case "CAMERA":
          case "BODY":
          case "SHADER":
            post += ":\n";
            break;
          case "TEMPLATE":
            post = "myTransform() {\n  return z => {\n    return {\n      ...z,\n      re: z.re,\n      im: z.im\n    }\n  }\n}\n";
            break;
          case "type":
            post += " n = ";
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
  name: "tutorial",
  bindKey: {
    win: "Ctrl-/|Ctrl-Shift-/",
    mac: "Command-/|Command-Shift-/"
  },
  exec: function() {
    tutorial();
  }
});
