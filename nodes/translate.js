const translate = {
  name: "translate",
  function: function(c){
    return function(z){
      let ans = compAdd(z, c);
      return {
        ...z,
        ...ans
      }
    }
  },
  parameters: [
    {
      name: "c",
      type: "complex"
    }
  ]
}
