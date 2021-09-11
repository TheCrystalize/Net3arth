const WIDTH = 800;
const HEIGHT = 800;

//initialize canvas
const c = document.getElementById("canvas");
const ctx = c.getContext("2d", {
  alpha: false
});
c.style.backgroundColor = "#12161b";
c.style.minWidth = WIDTH + 'px';
c.style.maxWidth = WIDTH + 'px';
c.width = WIDTH;
c.style.minHeight = HEIGHT + 'px';
c.style.maxHeight = HEIGHT + 'px';
c.height = HEIGHT;

ctx.imageSmoothingQuality = "high";

/**
 * Updates the canvas content
 */

function doStuff(){
  let val = {
    re: Math.random()-0.5,
    im: Math.random()-0.5,
    red: 1,
    green: 1,
    blue: 1
  };
  return val;
}

let buffer = Array(WIDTH * HEIGHT).fill([0,0,0]);
let img = new ImageData(WIDTH, HEIGHT);
for(let i=3;i<WIDTH*HEIGHT*4;i+=4){
  img.data[i]=255;
}

function draw() {
  let start = Date.now();
  do{
    for(let i = 0; i < 100000; i++){
      let val = doStuff();
      if(val.re > -0.5 && val.re < 0.5 && val.im > -0.5 && val.im < 0.5){
        let index = ((val.re+0.5)*WIDTH>>0) + ((-val.im+0.5)*WIDTH>>0)*HEIGHT;
        let t = buffer[index];
        buffer[index]=[t[0]+val.red,t[1]+val.green,t[2]+val.blue];
      }
    }
  } while(Date.now()-start < 12);

  for(let i=0;i<WIDTH*HEIGHT;i++){
    img.data[i*4] = buffer[i][0];
    img.data[i*4+1] = buffer[i][1];
    img.data[i*4+2] = buffer[i][2];
  }
  ctx.putImageData(img, 0, 0);

  //ctx.fillStyle = 'red';
  //ctx.fillRect(Math.random()*WIDTH, 0, WIDTH, HEIGHT);

  window.requestAnimationFrame(draw);
}

draw();
