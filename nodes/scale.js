const scale = {
  name: "scale",
  function: function(s){
    return function(z){
      let ans = compMultScalar(z, s);
      return {
        ...z,
        ...ans
      }
    }
  },
  parameters: [
    {
      name: "s",
      type: "float"
    }
  ]
}
