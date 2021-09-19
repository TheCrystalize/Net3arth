let tutorialState = 0;

let tutorials = [
`/*
 This is the Net3arth tutorial
 The goal of these tutorials are to teach you
 how to use Net3arth & program in 3arthLang.

 All of these tutorials are written in runnable
 3arthLang code. You can run the code by
 clicking the Run button in the top left,
 typing Ctrl-R (Cmd-R on Mac),
 or typing F1, searching for "run" and ENTER.

 Try it now!

 You should see a white circle appear.
 */

body: blurCircle()
        -> scale(0.5);
`,
`/*
  That was a fine circle, but it'll be
  much better once we start using some
  iteration. We do this most simply with "choose".
*/

body:
choose{
  1: blurCircle()
  -> circleInv();
  1: rotate(0.5)
  -> arcsinh()
  -> scale(0.7853981633974)
  -> trigTanh();
};
`,

];

function tutorial() {
  let pos = editor.getCursorPosition();

  editor.setValue(tutorials[tutorialState]);
  editor.moveCursorToPosition({row:0, pos:0});
  editor.clearSelection();

  tutorialState = (tutorialState+1)%tutorials.length;
}
