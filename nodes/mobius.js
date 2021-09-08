const mobius = {
  name: "mobius",
  function: function(mobius){
    return function(z, a, b, c, d){
      z = (a * z + c) / (b * z + d);
      return z;
    }
  }
}
