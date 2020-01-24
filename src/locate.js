import { minimize_Powell } from 'optimization-js';
import minimize_cobyla from './3rdparty/jscobyla/index';

import { E, Point, Circle } from './geometry';

const sum = vals => vals.reduce((agg, v) => agg + v, 0);

const Norm = (x, y, mode='2D') => {
  if (mode === '2D') return ((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2) ** .5;
  if (mode === '3D') return ((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2) ** .5;
  if (mode === 'Earth1') return E.gcd(x[0], x[1], y[0], y[1]);
};

const sumError = (x, c, r, mode) => {
  const l = c.length;
  let e = 0;
  for (let i = 0; i < l; i++) {
    e += (Norm(x, c[i].std(), mode) - r[i]) ** 2;
  }
  return e;
};

const isDisjoint = (cA, fg = 0) => {
  // Currently this function does not support sophisticated checking
  // for  disjoint area
  // on earth surface models
  const l = cA.length;
  for (let i = 0; i < l; i++) {
    for (let j = i+1; j < l; j++) {
      if (!cA[j].touch(cA[i], fg)) return true;
    }
  }
  return false;
};

function lse(cA, mode='2D', cons = true) {
  const l = cA.length;
  const r = cA.map(w => w.r);
  const c = cA.map(w => w.c);
  const S = sum(r);
  const W = r.map(w => (S - w) / ((l - 1) * S));

  let p0 = new Point(0, 0, 0); // Initialized point
  for (let i = 0; i < l; i++) {
    p0 = p0.sum(W[i]).sum(c[i]);
  }

  const x0 = (mode === '2D' || mode === 'Earth1')
    ? [p0.x, p0.y]
    : [p0.x, p0.y.p0.z];

  const fg1 = mode === 'Earth1' ? 1 : 0;

  let ans;
  if (cons) {
    console.info('GC-LSE geolocating...');
    if (!isDisjoint(cA, fg1)) {
      const constrain = (x, beaconIdx) => r[beaconIdx] - Norm(x, c[beaconIdx].std(), mode);

      let x = x0.slice(); // result

      const res = minimize_cobyla(
        (n, m, x, con) => { // objective function
          for (let i = 0; i < m; i++) {
            con[i] = constrain(x, i); // one constrain per beacon
          }
          return sumError(x, c, r, mode);
        },
        x.length, // # vars
        l, // # constraints
        x0.slice(), // result
        1, // rhobeg
        1e-5, // rhoend
        0, // iprint
        1000 // maxfun
      );

      console.log({ res, x });

      ans = x;

      // const res = fmin_cobyla(sum_error, x0, cL, args=(c, r, mode), consargs=(), rhoend=1e-5)
    } else {
      throw new Error('Disjoint');
    }
  } else {
    console.info('LSE geolocating...');

    const res = minimize_Powell(
      x => sumError(x, c, r, mode),
      x0
    );

    ans = res.argument;
  }

  return new Point(...ans);
}

function locate(beacons, { mode = '2D', constrain = true }) {
  if (mode !== '2D' && mode !== '3D' && mode !== 'Earth1') {
    throw new Error(`Mode not supported: ${mode}`);
  }

  const circles = beacons.map(({ distance, ...coords }) => {
    const c = mode === '2D'
      ? new Point(coords.x, coords.y)
      : mode === '3D'
        ? new Point(coords.x, coords.y, coords.z)
        : new Point(coords.lng, coords.lat);

    return new Circle(c, distance); // in Earth mode, gcd distances are specified in meters
  });

  const { x, y, z } = lse(circles, mode, constrain);

  if (mode === 'Earth1') {
    let lng = x;
    let lat = y;

    // make lat, lngs uniform
    while (lng > 180) { lng -= 360 }
    while (lng < -180) { lng += 360 }
    while (lat > 90) { lat -= 180 }
    while (lat < -90) { lat += 180 }

    return { lat, lng };
  }

  return mode === '2D'
    ? { x, y }
    : { x, y, z };
}

export default locate;
