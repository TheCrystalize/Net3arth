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

function sub(z, c) {
  return {
    re: z.re - c.re,
    im: z.im - c.im
  }
}

function addScalar(z, s) {
  return {
    re: z.re + s,
    im: z.im
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

function sqrt(z) {
  const s = Math.sqrt(z.re * z.re + z.im * z.im),
    sgn = z.im < 0 ? -1 : 1;
  return multScalar({
    re: Math.sqrt(s + z.re),
    im: sgn * Math.sqrt(s - z.re)
  }, 0.5 * Math.SQRT2 );
}

function log(z) {
  return {
    re: 0.5 * Math.log(z.re * z.re + z.im * z.im),
    im: Math.atan2(z.im, z.re)
  }
}

function pow(z, p) {
  const n = p * Math.atan2(z.im, z.re);
  return {
    re: Math.cos(n),
    im: Math.sin(n) * Math.exp(p * Math.log(z.re * z.re + z.im * z.im))
  }
}

function exp(z) {
  const e = Math.exp(z.re);
  return {
    re: Math.cos(z.im) * e,
    im: Math.sin(z.im) * e
  }
}

/*Transforms*/
function arcsinh(){
  return function(z){
    let ans = multScalar(
      log(
        add(
          z,
          sqrt(
            addScalar(mult(z, z), 1)
          )
        )
      ),
      2 / Math.PI);
    return {
      ...z,
      ...ans
    }
  }
}

function splits(x, y) {
  return function(z) {
    let xoff = z.re >= 0 ? x : -x,
      yoff = z.im >= 0 ? y : -y;
    let ans = {
      re: z.re + xoff,
      im: z.im + yoff
    }
    return {
      ...z,
      ...ans
    }
  }
}

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

function translate(real, imaginary) {
  return function(z) {
    return {
      ...z,
      re: z.re + real,
      im: z.im + imaginary
    }
  }
}

function rotate(theta) {
  return function(z) {
    return {
      ...z,
      re: - Math.cos(theta) * z.re - Math.sin(theta) * z.im,
      im: - Math.sin(theta) * z.re + Math.cos(theta) * z.im
    }
  }
}

function blurCircle(){
  return function(z) {
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
}

const BUILT_IN_TRANSFORMS = {
  arcsinh: arcsinh,
  splits: splits,
  mobius: mobius,
  scale: scale,
  translate: translate,
  rotate: rotate,
  blurCircle: blurCircle,
};

const BUILT_IN_TRANSFORMS_PARAMS = {
  arcsinh: [],
  splits: [
    {
      name: "x",
      type: "number",
      default: 0
    },
    {
      name: "y",
      type: "number",
      default: 0
    }
  ],
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
  translate: [
    {
      name: "real",
      type: "number",
      default: 0
    },
    {
      name: "imaginary",
      type: "number",
      default: 0
    }
  ],
  rotate: [
    {
      name: "theta",
      type: "number",
      default: 0
    }
  ],
  blurCircle: []
};
