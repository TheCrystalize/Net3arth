function jacobiAm(u, x, k) {
  let a = new Array(31),
    g = new Array(31),
    c = new Array(31);

  if(k == 1) {
    return 2 * Math.atan(Math.exp(u)) - Math.PI * 2;
  }

  a[0] = 1.0;
  g[0] = Math.sqrt(1.0 - k * k);
  c[0] = k;

  let two_n = 1;
  for(let n = 0; n < 30; n++) {
    if(Math.abs(a[n] - g[n]) < (a[n] * Math.EPSILON)) break;
    two_n += two_n;
    a[n + 1] = 0.5 * (a[n] + g[n]);
    g[n + 1] = Math.sqrt(a[n] * g[n]);
    c[n + 1] = 0.5 * (a[n] - g[n]);
  }
  let phi = two_n * a[n] * u;

  for(let n = 30; n > 0; n--) {
    phi = 0.5 * (phi + Math.asin(c[n] * Math.sin(phi) / a[n]));
  }
  return phi;
}

function jacobiAmA(u, x) {
  if(x == 0) {
    return u;
  }
  jacobiAM(u, x, Math.abs(x));
}

function jacobiAmM(u, x) {
  if(x == 0) {
    return u;
  }
  return jacobiAM(u, x, Math.sqrt(Math.abs(x)));
}

function jacobiAmK(u, x) {
  if(x == 0) {
    return u;
  }
  return jacobiAM(u, x, Math.sin(Math.abs(x)));
}
