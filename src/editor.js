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
      ["keyword", ["choose", "xaos"]],
      ["type", ["complex", "number", "string", "bool", "array", "object"]],
      ["transform", StandardLibNamesArray],
      ["constructor", StandardLibConstructorNamesArray],
      ["helper", StandardLibHelperNamesArray]
    ];
    callback(null, wordList.map(function(words) {
      return words[1].map(function(word) {
        let post = word;
        switch (words[0]) {
          case "CAMERA":
          case "BODY":
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
            if(BUILT_IN_TRANSFORMS_PARAMS.hasOwnProperty(word)) {
              post += "(";
              for(let i = 0; i < BUILT_IN_TRANSFORMS_PARAMS[word].length; i++) {
                switch (BUILT_IN_TRANSFORMS_PARAMS[word][i].type) {
                  case "number":
                  case "bool":
                    post += BUILT_IN_TRANSFORMS_PARAMS[word][i].default;
                    break;
                  case "array":
                  case "object":
                    post += JSON.stringify(BUILT_IN_TRANSFORMS_PARAMS[word][i].default);
                    break;
                  case "string":
                    post += '"' + BUILT_IN_TRANSFORMS_PARAMS[word][i].default+'"';
                    break;
                  case "complex":
                    if(BUILT_IN_TRANSFORMS_PARAMS[word][i].default.im === 0) {
                      post += BUILT_IN_TRANSFORMS_PARAMS[word][i].default.re;
                    } else if(BUILT_IN_TRANSFORMS_PARAMS[word][i].default.re === 0) {
                      post += BUILT_IN_TRANSFORMS_PARAMS[word][i].default.im + 'i';
                    } else if(BUILT_IN_TRANSFORMS_PARAMS[word][i].default.im < 0) {
                      post += BUILT_IN_TRANSFORMS_PARAMS[word][i].default.re + "-" + (-BUILT_IN_TRANSFORMS_PARAMS[word][i].default.im) + 'i';
                    } else {
                      post += BUILT_IN_TRANSFORMS_PARAMS[word][i].default.re + "+" + BUILT_IN_TRANSFORMS_PARAMS[word][i].default.im + 'i';
                    }
                    break;
                }
                if(i < BUILT_IN_TRANSFORMS_PARAMS[word].length - 1) {
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
    win: "Ctrl-R",
    mac: "Command-R"
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
