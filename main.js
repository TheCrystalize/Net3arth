const WIDTH = 800;
const HEIGHT = 800;

//initialize canvas
const c = document.getElementById("canvas");
const ctx = c.getContext("2d", {
  alpha: false
});
c.style.backgroundColor = "red";
c.width = WIDTH;
c.height = HEIGHT;

/**
 * Updates the canvas content
 */
function draw() {
  ctx.imageSmoothingQuality = "high";
  ctx.fillStyle = "#224488";
  ctx.fillRect(0, 0, WIDTH, HEIGHT);
  ctx.fillStyle = "black";
  ctx.fillRect(WIDTH/4, HEIGHT/4, WIDTH/2, HEIGHT/2);
}

draw();
