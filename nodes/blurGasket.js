const blurGasket = {
  name: "blurGasket",
  function: function(){
    return function(z){
      let a = Math.random() * Math.PI * 2,
          r = 1 / Math.sqrt(Math.random() - Math.random());
      let ans = {
        re: Math.random() - 0.5,
        im: Math.sin(a) * r
      }
      return {
        ...z,
        ...ans
      }
    }
  },
}
