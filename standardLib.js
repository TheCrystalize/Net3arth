/*helper functions*/
function C(real, imaginary) {
  return {re: real, im: imaginary};
}

function div(z, c) {
  const s = 1 / (c.re * c.re + c.im * c.im);
  return {
    re: (z.re * c.re + z.im * c.im) * s,
    im: (z.im * c.re - z.re * c.im) * s
  }
}

function add(z, c) {
  return {
    re: z.re + c.re,
    im: z.im + c.im
  }
}

function mult(z, c) {
  return {
    re: z.re * c.re - z.im * c.im,
    im: z.re * c.im + z.im * c.re
  }
}

function multScalar(z, s) {
  return {
    re: z.re * s,
    im: z.im * s
  }
}
/*Transforms*/
function mobius(a, b, c, d) {
  return function(z) {
    let ans = div(add(mult(a, z), b), add(mult(c, z), d));
    return {
      ...z,
      ...ans
    }
  }
}

function scale(s) {
  return function(z) {
    let ans = multScalar(z, s);
    return {
      ...z,
      ...ans
    }
  }
}

function blurCircle(z){
  let a = Math.random() * Math.PI * 2,
      r = Math.sqrt(Math.random());
  let ans = {
    re: Math.cos(a) * r,
    im: Math.sin(a) * r
  }
  return {
    ...z,
    ...ans
  }
}

const BUILT_IN_TRANSFORMS = {
  mobius: [
    {
      name: "a",
      type: "complex",
      default: {re:1,im:0}
    },
    {
      name: "b",
      type: "complex",
      default: {re:0,im:0}
    },
    {
      name: "c",
      type: "complex",
      default: {re:0,im:0}
    },
    {
      name: "d",
      type: "complex",
      default: {re:1,im:0}
    }
  ],
  scale: [
    {
      name: "s",
      type: "number",
      default: 1
    }
  ],
  blurCircle: []
};
