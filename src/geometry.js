const res = 0.0001; // checking accuracy

class Point {
  constructor(...args) {
    const l = args.length;

    this.dim = l === 1 ? args[0].length : l;
    if (l === 1) {
      this.x = args[0][0];
      this.y = args[0][1];
      this.z = this.dim > 2 ? args[0][2] : 0;
    } else {
      this.x = args[0];
      this.y = args[1];
      this.z = this.dim > 2 ? args[2] : 0;
    }
  }

  sum(other) {
    if (other instanceof Point) {
      return new Point(this.x + other.x, this.y + other.y, this.z + other.z);
    } else {
      return new Point(this.x + other, this.y + other, this.z + other);
    }
  }

  sub(other) {
    if (other instanceof Point) {
      return new Point(this.x - other.x, this.y - other.y, this.z - other.z);
    } else {
      return new Point(this.x - other, this.y - other, this.z - other);
    }
  }

  mult(other) {
    if (other instanceof Point) {
      return new Point(this.x * other.x, this.y * other.y, this.z * other.z);
    } else {
      return new Point(this.x * other, this.y * other, this.z * other);
    }
  }

  div(other) {
    if (other instanceof Point) {
      return new Point(this.x / other.x, this.y / other.y, this.z / other.z);
    } else {
      return new Point(this.x / other, this.y / other, this.z / other);
    }
  }

  dist(other) {
    return ((this.x - other.x) ** 2 + (this.y - other.y) ** 2 + (this.z - other.z) ** 2) ** 0.5;
  }

  std() {
    return this.dim === 2
      ? [this.x, this.y]
      : [this.x, this.y, this.z];
  }
}

class Circle {
  constructor(p, r) {
    this.c = p; // Point
    this.r = +r; // num
  }

  touch(o, earthSurface = false) {
    const d = this.c.dist(o.c);
    const met = !earthSurface
      ? this.r + o.r
      : 180 * (this.r + o.r) / (Earth.R * Math.PI);

    return d < (met + res);
  }
}

class Sphere {
  constructor(c, R) {
    this.c = c;
    this.R = R;
  }

  gcd(lon1, lat1, lon2, lat2) {
    // Calculate the great circle distance (in meters) between two points on the earth (specified in lng/lat degrees)

    // convert to radians
    const [_lon1, _lat1, _lon2, _lat2] = [lon1, lat1, lon2, lat2].map(a => a * Math.PI / 180);

    // haversine formula
    const dlon = _lon2 - _lon1;
    const dlat = _lat2 - _lat1;
    const a = Math.sin(dlat / 2) ** 2 + Math.cos(_lat1) * Math.cos(_lat2) * Math.sin(dlon / 2) ** 2;
    const c = 2 * Math.asin(Math.sqrt(a));

    return Earth.R * c;
  }
}

const Earth = new Sphere(new Point(0, 0, 0), 6378100); // radius in meters

export { Point, Circle, Earth };