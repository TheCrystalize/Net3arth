let tutorialState = -1;

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

body:
blurCircle()
-> scale(0.5);

/*
 The above code is defining the body of the code
 and applying a circle function which is piped
 into a scale function using the "->" operator.

 ~continue to the next page~
   ~Press Ctrl-? or Cmd-?~
 */
`,
  `/*
 In addition to the body, there's the camera
 The camera is run before the points are
 drawn, and then discarded.
 */

body:
blurCircle();

camera:
scale(0.5);

/*
 ~continue to the next page~
   ~Press Ctrl-? or Cmd-?~
 */
`,
  `/*
 we can define a new transform like this:
transform <name>(<type> <name>,...):
  <transforms>;

 Here's an example that moves points
 halfway to a given coordinate:
 */
transform halfway(number x, number y):
  scale(1/2) -> translate(x/2, y/2);

/*
 this is the main (if not only) transform
 used in most introductions to the chaos game;
 https://www.google.com/search?q=chaos+game

 And the chaos game forms the foundation on
 which this all works.
 */

/*
 As a simple demonstration, here's a circle
 that starts at (-0.4,-0.4) and moves towards (0.49,0.5)
 */
body:
choose{
  1: halfway(0.49, 0.5);

  // this just draws a circle in the top left
  1: blurCircle()
    -> scale(1/20)
    -> translate(-0.4, -0.4);
};

/*
 this is a shader. It runs on the pixel level
 and lets us do things with the rendered image.
 In this case, we're just doing some gamma
 correction so that the circles are brighter
 */
shader:
gamma(10);
`,
`transform halfway(number x, number y):
  scale(1/2) -> translate(x/2, y/2);

/*
 if we add another halfway to our choose,
 choose will randomly pick between them
 and the result is a line.
 */

body:
choose{
  1: halfway(0.49, 0.5);
  1: halfway(0.49, -0.5);

  // this just draws a circle in the top left
  1/4:
    blurCircle()
    -> scale(1/20)
    -> translate(-0.4, 0);
};
`,
`transform halfway(number x, number y):
  scale(1/2) -> translate(x/2, y/2);

/*
 if we add one more halfway to our choose,
 choose will randomly we will get lines
 connecting all 3 of our points;
 creating a triangle, but also a patern in the middle,
 which I find most easily understood via this animation:
 https://upload.wikimedia.org/wikipedia/commons/b/b6/Sierpinski_Chaos.gif

 I added a color() function to match that animations' coloring.
 It sets the color to the RGB values specified.
 */

body:
choose{
  1: halfway(0.49, 0.5)  -> color(colorRGB(0, .5, 1));
  1: halfway(0.49, -0.5) -> color(colorRGB(1,  0, 0));
  1: halfway(-0.4, 0)    -> color(colorRGB(1,  1, 0));

  // this just draws a circle in the top left
  1:
    reset() // resets coordinates, color, and alpha
    -> blurCircle()
    -> scale(1/20)
    -> translate(-0.4, 0);
};

shader:
gamma(3);
// sometimes higher gammas make details
// easier to see, but it adds graininess
`,
`
/*
 You can use Ctrl+Shift+Space (Cmd+Shift+Space for Mac)
 to see a list of all the functions.
 F1 to see a list of keyboard shortcuts.
 Ctrl-comma to open the editor settings.
 */
`
];

function tutorial(back = false) {
  if(tutorialState < 0){
    tutorialState = 0;
  }
  else{
    tutorialState = (tutorials.length + tutorialState + (back ? -1 : 1)) % tutorials.length;
  }

  let pos = editor.getCursorPosition();

  editor.setValue(tutorials[tutorialState]);
  editor.moveCursorToPosition({
    row: 0,
    pos: 0
  });
  editor.clearSelection();
}
