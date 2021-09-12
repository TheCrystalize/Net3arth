function blurCircle(){
  return function(z){
    let a = Math.random() * Math.PI * 2,
        r = Math.sqrt(Math.random());
    let ans = {
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    };
    return {
      ...z,
      ...ans
    };
  }
}

const blurCircleData = {
  name: "blurCircle",
  parameters: []
}
