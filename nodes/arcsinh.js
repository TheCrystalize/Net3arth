function arcsinh(){
  return function(z){
    let ans = (2 / Math.PI) * log(z + sqrt(z * z + 1.0));
    return {
      ...z,
      ...ans
    };
  }
}

const arcsinhData = {
  name: "arcsinh",
  parameters: []
}
