'use strict';

/**
 * @namespace
 * @type {Object} Visualization object
 */
var vz = {
  name      : 'CFD is Colorful Fluid Dynamics',
  type      : 'visualization',
  tags      : ['CFD', 'canvas'],
  screen    : null,
  canvas    : null,
  context   : null,
  width     : 0,
  height    : 0,
  bands     : [],
  initialized: false,
  // begin application specific code
  particles : 30000,
  px        : null, // particle x locations
  py        : null, // particle y locations
  pc        : null, // particle hue [0, 1]
  pl        : null, // particle life
  lifetime  : 100, // particle lifetime
  field     : null, // FluidField
  fieldDisplay: null, // FluidFieldDisplay
  detector  : null, // BeatDetector
  showVelocity: false, // toggle to show/hide velocity field
  showParticles: true, // toggle to show/hide particles
  radius    : 8, // size of particle spawning ring
  startTime : 0, // start time
  now       : 0, // current time
  frames    : 0, // number of frames rendered
  options: { },
  audio: {
    audio: function (event) {
      vz.detector.process(event.audio);
      return true;
    },
    pause: function (event) {
      return true;
    },
    reset: function (event) {
      return true;
    }
  }
};

/** 
 * Reset particle position, e.g. after it dies
 * @param i the index of particle to be reset
 */
vz.spawn = function(i) {
  var t = i / vz.particles;
  var r = vz.radius + Math.random();
  vz.px[i] = vz.field.width * 0.5 + r * Math.cos(t * 2 * Math.PI);
  vz.py[i] = vz.field.height * 0.5 + r * Math.sin(t * 2 * Math.PI);
  vz.pc[i] = t;
  vz.pl[i] = vz.lifetime;
};

/**
 * Resets the lifetime of particles of a certain color, triggered on strong beats
 * @param min the minimum hue [0, 1]
 * @param max the maximum hue [0, 1]
 */
vz.pulse = function(min, max) {
  for (var i = 0; i < vz.particles; i++) {
    if (Math.random() < 0.5 && min < vz.pc[i] && vz.pc[i] < max) {
      // subtract random number so particles don't all die at once
      vz.pl[i] = vz.lifetime - Math.floor(Math.random() * 15);
    }
  }
};

/**
 * Animate! 
 */
vz.animate = function() {
  
  if (!vz.initialized) return;
  requestAnimationFrame(vz.animate);

  vz.now = new Date();
  vz.time = vz.now - vz.startTime;
  
  vz.updateVelocities(); 
  vz.field.update();

  // update particle positions and lifetime
  var jitter, vx, vy;
  for (var i = 0; i < vz.particles; i++) {
    // add jitter to old particles for diffuse smoke effect
    jitter = (1 - vz.pl[i] / vz.lifetime);
    vx = 10 * vz.field.getXVelocity(vz.px[i], vz.py[i]);
    vy = 10 * vz.field.getYVelocity(vz.px[i], vz.py[i]);
    vz.px[i] += vx + (Math.random() - 0.5) * jitter;
    vz.py[i] += vy + (Math.random() - 0.5) * jitter;
    vz.pl[i]--;
    if (vz.pl[i] < 1 || vz.px[i] < 1 || vz.px[i] > vz.field.width || vz.py[i] < 1 || vz.py[i] > vz.field.height - 1) {
      vz.spawn(i);
    }
  }

  if (vz.showParticles)
    vz.display.renderParticles(vz.field, vz.px, vz.py, vz.pc, vz.pl);
  
  if (vz.showVelocity)
    vz.display.renderVelocityField(vz.field);

};

vz.updateVelocities = function() {
  var explode = 1.5;
  var jet = 1.2;
  var bigJet = 1.5;
  var maxVelocity = 4;
  console.log("band likelihood: " + vz.detector.bandLikelihood);
  var offset = vz.now * 0.002;
  var x, y, theta;
  vz.radius = (vz.detector.beatLikelihood > explode) ? 15 : 8;
  for (var i = 0; i < vz.detector.bands; i++) {
    theta = i / vz.detector.bands * 2 * Math.PI + offset;
    x = Math.floor(vz.field.width / 2 + vz.radius * Math.cos(theta));
    y = Math.floor(vz.field.height / 2 + vz.radius * Math.sin(theta));
    var likelihood = vz.detector.bandLikelihood[i];
    if (likelihood > jet) {
      var v = likelihood * 2;
      if (v > maxVelocity) v = maxVelocity;
      var vx = vz.field.getXVelocity(x, y) + v * Math.cos(theta);
      var vy = vz.field.getYVelocity(x, y) + v * Math.sin(theta);
      vz.field.setVelocity(x, y, vx, vy);
      if (bigJet < likelihood) {
        var min = (i / vz.detector.bands) - 0.01;
        var max = min + 0.02;
        vz.pulse(min, max);
      }
    } 
  }  
};

/**
 * Start visualization
 * @param options
 */
vz.start = function (options) {

  vz.screen = options.screen;
  vz.ww = options.width;
  vz.wh = options.height;

  vz.canvas = document.createElement('canvas');
  vz.screen.appendChild(vz.canvas);

  vz.canvas.width = options.width - 20;
  vz.canvas.height = options.height - 20; 

  vz.context = vz.canvas.getContext('2d');
  vz.context.globalCompositeOperation = 'destination-over';

  vz.px = new Float32Array(vz.particles);
  vz.py = new Float32Array(vz.particles);
  vz.pc = new Float32Array(vz.particles);
  vz.pl = new Int16Array(vz.particles);

  vz.field = new FluidField(vz.canvas);
  vz.display = new FluidFieldDisplay(vz.field, vz.canvas);
  vz.detector = new BeatDetector();

  vz.startTime = new Date();

  for (var i = 0; i < vz.particles; i++) {
    vz.spawn(i);
    vz.pl[i] = Math.floor(Math.random() * vz.lifetime);
  }

  vz.initialized = true;

  requestAnimationFrame(vz.animate);
};

/**
 * Stop visualization
 */
vz.stop = function () {
  if (!vz.initialized) return;
  vz.screen.removeChild(vz.canvas);
  vz.canvas = null;
  vz.context = null;
};

/**
 * Start visualization 
 * @param options 
 */
vz.fadeIn = function (options, step) {
  vz.start(options);
  vz.context.globalAlpha = 1;
};

/** Stop visualization fading out */
vz.fadeOut = function (step) {
  vz.context.globalAlpha = 0;
  vz.stop();
};

/**
 * Resize the visualization
 * @param width
 * @param height
 */
vz.resize = function (width, height) {
  if (!vz.initialized) return;
  vz.canvas.width = width - 20;
  vz.canvas.height = height - 20;
};

// Export API
exports.name    = vz.name;
exports.type    = vz.type;
exports.tags    = vz.tags;
exports.start   = vz.start;
exports.stop    = vz.stop;
exports.fadeIn  = vz.fadeIn;
exports.fadeOut = vz.fadeOut;
exports.resize  = vz.resize;
exports.audio   = vz.audio;


function BeatDetector() {

  this.integration = 10;

  this.bands = 10;
  this.bandHistory = [];
  this.bandLikelihood = [];

  this.levelHistory = [];
  this.beatLikelihood = 0;

  for (var i = 0; i < this.bands; i++) {
    this.bandHistory[i] = []; 
    this.bandLikelihood[i] = 0;
  }

}


BeatDetector.prototype.process = function(audio) {
  
  // compute sample mean
  function computeMean(arr) {
    var sum = 0;
    for (var i = 0; i < arr.length; i++) {
      sum += arr[i];
    }
    return sum / arr.length;
  }

  // compute sample standard deviation
  function computeStdev(arr, mean) {
    var sum = 0;
    for (var i = 0; i < arr.length; i++) {
      sum += (arr[i] - mean) * (arr[i] - mean);
    }
    return Math.sqrt(sum / arr.length);
  }

  // convert audio.spectrum into mono for easier processing
  var spectrum = [];
  for (var i = 0; i < audio.spectrum.left.length; i++) {
    spectrum[i] = (parseFloat(audio.spectrum.left[i]) + parseFloat(audio.spectrum.right[i])) / 2;
    // spectrum[i] = Math.exp(spectrum[i]);
  }

  // get average level weighted to emphasize low freq
  var average = 0;
  for (var i = 0; i < spectrum.length; i++) {
    var weight = 2 - (i / spectrum.length);
    average += weight * spectrum[i];
  }
  average /= spectrum.length;

  // if we have fully populated integration window, do analysis
  if (this.levelHistory.length == this.integration) {
    var mean = computeMean(this.levelHistory);
    var stdev = computeStdev(this.levelHistory, mean);
    this.beatLikelihood = Math.abs(average - mean) / (stdev + 0.01);
    this.levelHistory.shift();
  } 
  // add current average level to history
  this.levelHistory.push(average);

  // analyze per band
  var bandsPerBin = spectrum.length / this.bands; // number of spectrum bands per analysis band
  for (var i = 0; i < this.bands; i++) {
    var current = 0;
    for (var j = 0; j < bandsPerBin; j++) {
      current += spectrum[i * bandsPerBin + j] / bandsPerBin;
    }
    // if we have a full integration window, do analysis
    if (this.bandHistory[i].length == this.integration) {
      var mean = computeMean(this.bandHistory[i]);
      var stdev = computeStdev(this.bandHistory[i], mean);    
      console.log(mean + " " + stdev);

      this.bandLikelihood[i] = Math.abs(current - mean) / (stdev + 0.01);
      this.bandHistory[i].shift();
    }
    // add current band level to history
    this.bandHistory[i].push(current);
    console.log(this.bandHistory[i].length);
  }
  
}

function FluidField() {
  this.iterations = 10;   // solver iterations
  this.dt = 0.1;          // timestep
  this.u = null;          // velocity x field
  this.u_prev = null;
  this.v = null;          // velocity y field
  this.v_prev = null;
  this.width = 120;       // mesh size x
  this.height = 80;       // mesh size y
  this.rowSize;           // rowstride
  this.size;              // mesh points
  this.interpolate = true;// interpolate when setting/getting fields
  this.damp = 0.99;       // velocity damping
  this.reset();
}

FluidField.prototype.reset = function() {
  this.rowSize = this.width + 2;
  this.size = (this.width + 2) * (this.height + 2);
  this.dens = new Float32Array(this.size);
  this.dens_prev = new Float32Array(this.size); 
  this.u = new Float32Array(this.size);
  this.u_prev = new Float32Array(this.size); 
  this.v = new Float32Array(this.size);
  this.v_prev = new Float32Array(this.size);
}

FluidField.prototype.addFields = function (x, s, dt) {
  for (var i = 0; i < this.size; i++) {
    x[i] += s[i] * dt;
  }
}

FluidField.prototype.enforceBCs = function(b, x) {
  var width = this.width;
  var height = this.height;
  var rowSize = this.rowSize;
  if (b == 1) {
      for (var i = 1; i <= width; i++) {
          x[i] =  x[i + rowSize];
          x[i + (height+1) * rowSize] = x[i + height * rowSize];
      }
      for (var j = 1; i <= height; i++) {
          x[j * rowSize] = -x[1 + j * rowSize];
          x[(width + 1) + j * rowSize] = -x[width + j * rowSize];
      }
  } else if (b == 2) {
      for (var i = 1; i <= width; i++) {
          x[i] = -x[i + rowSize];
          x[i + (height + 1) * rowSize] = -x[i + height * rowSize];
      }
      for (var j = 1; j <= height; j++) {
          x[j * rowSize] =  x[1 + j * rowSize];
          x[(width + 1) + j * rowSize] =  x[width + j * rowSize];
      }
  } else {
      for (var i = 1; i <= width; i++) {
          x[i] =  x[i + rowSize];
          x[i + (height + 1) * rowSize] = x[i + height * rowSize];
      }
      for (var j = 1; j <= height; j++) {
          x[j * rowSize] =  x[1 + j * rowSize];
          x[(width + 1) + j * rowSize] =  x[width + j * rowSize];
      }
  }
    var maxEdge = (height + 1) * rowSize;
    x[0]                      = 0.5 * (x[1] + x[rowSize]);
    x[maxEdge]                = 0.5 * (x[1 + maxEdge] + x[height * rowSize]);
    x[(width+1)]         = 0.5 * (x[width] + x[(width + 1) + rowSize]);
    x[(width+1)+maxEdge] = 0.5 * (x[width + maxEdge] + x[(width + 1) + height * rowSize]);
}

FluidField.prototype.linSolve = function(b, x, x0, a, c) {
  var width = this.width;
  var height = this.height;
  var rowSize = this.rowSize;
  if (a == 0 && c == 1) {
    for (var j = 1; j <= height; j++) {
      var row = j * rowSize; 
      ++row;
      for (var i = 0; i < width; i++) {
        x[row] = x0[row];
        ++row;
      }
    }
    this.enforceBCs(b, x);
  } else { 
    var cinv = 1 / c;
    for (var k = 0; k < this.iterations; k++) {
      for (var j = 1; j <= height; j++) {
        var last = (j - 1) * rowSize;
        var current = j * rowSize;
        var next = (j + 1) * rowSize;
        var lastx = x[current];
        ++current;
        for (var i = 1; i < width; i++) {
          lastx = x[current] = (x0[current] + a * (lastx + x[++current] + x[++last] + x[++next])) * cinv;
        }
      }
      this.enforceBCs(b, x);
    }
  }
}

FluidField.prototype.diffuse = function(b, x, x0) {
  this.linSolve(b, x, x0, 0, 1);
}

FluidField.prototype.linSolve2 = function(x, x0, y, y0, a, c) {
    if (a == 0 && c == 1) {
        for (var j = 1 ; j <= this.height; j++) {
            var current = j * this.rowSize;
            ++current;
            for (var i = 0; i < this.width; i++) {
                x[current] = x0[current];
                y[current] = y0[current];
                ++current;
            }
        }
        this.enforceBCs(1, x);
        this.enforceBCs(2, y);
    } else {
        var invC = 1 / c;
        for (var k = 0 ; k < this.iterations; k++) {
            for (var j=1 ; j <= this.height; j++) {
                var last = (j - 1) * this.rowSize;
                var current = j * this.rowSize;
                var next = (j + 1) * this.rowSize;
                var lastX = x[current];
                var lastY = y[current];
                ++current;
                for (var i = 1; i <= this.width; i++) {
                    lastX = x[current] = (x0[current] + a * (lastX + x[current] + x[last] + x[next])) * invC;
                    lastY = y[current] = (y0[current] + a * (lastY + y[++current] + y[++last] + y[++next])) * invC;
                }
            }
            this.enforceBCs(1, x);
            this.enforceBCs(2, y);
        }
    }
}

FluidField.prototype.diffuse2 = function(x, x0, y, y0, dt) {
  this.linSolve2(x, x0, y, y0, 0, 1);
}

FluidField.prototype.advect = function(b, d, d0, u, v, dt) {
  var rowSize = this.rowSize;
  var Wdt0 = this.dt * this.width;
  var Hdt0 = this.dt * this.height;
  var Wp5 = this.width + 0.5;
  var Hp5 = this.height + 0.5;
  for (var j = 1; j<= this.height; j++) {
      var pos = j * rowSize;
      for (var i = 1; i <= this.width; i++) {
          var x = i - Wdt0 * u[++pos]; 
          var y = j - Hdt0 * v[pos];
          if (x < 0.5)
              x = 0.5;
          else if (x > Wp5)
              x = Wp5;
          var i0 = x | 0;
          var i1 = i0 + 1;
          if (y < 0.5)
              y = 0.5;
          else if (y > Hp5)
              y = Hp5;
          var j0 = y | 0;
          var j1 = j0 + 1;
          var s1 = x - i0;
          var s0 = 1 - s1;
          var t1 = y - j0;
          var t0 = 1 - t1;
          var row1 = j0 * rowSize;
          var row2 = j1 * rowSize;
          d[pos] = s0 * (t0 * d0[i0 + row1] + t1 * d0[i0 + row2]) + s1 * (t0 * d0[i1 + row1] + t1 * d0[i1 + row2]);
      }
  }
  this.enforceBCs(b, d);
}

FluidField.prototype.project = function(u, v, p, div) {
  var width = this.width;
  var height = this.height;
  var rowSize = this.rowSize;
  var h = -0.5 / Math.sqrt(width * height);
  for (var j = 1 ; j <= height; j++ ) {
      var row = j * rowSize;
      var previousRow = (j - 1) * rowSize;
      var prevValue = row - 1;
      var currentRow = row;
      var nextValue = row + 1;
      var nextRow = (j + 1) * rowSize;
      for (var i = 1; i <= width; i++ ) {
          div[++currentRow] = h * (u[++nextValue] - u[++prevValue] + v[++nextRow] - v[++previousRow]);
          p[currentRow] = 0;
      }
  }
  this.enforceBCs(0, div);
  this.enforceBCs(0, p);
  
  this.linSolve(0, p, div, 1, 4 );
  var wScale = 0.5 * width;
  var hScale = 0.5 * height;
  for (var j = 1; j<= height; j++ ) {
      var prevPos = j * rowSize - 1;
      var currentPos = j * rowSize;
      var nextPos = j * rowSize + 1;
      var prevRow = (j - 1) * rowSize;
      var currentRow = j * rowSize;
      var nextRow = (j + 1) * rowSize;
      for (var i = 1; i<= width; i++) {
          u[++currentPos] -= wScale * (p[++nextPos] - p[++prevPos]);
          v[currentPos]   -= hScale * (p[++nextRow] - p[++prevRow]);
      }
  }
  this.enforceBCs(1, u);
  this.enforceBCs(2, v);
}

FluidField.prototype.velocityStep = function(u, v, u0, v0, dt) {
    this.addFields(u, u0, dt);
    this.addFields(v, v0, dt);
    var temp = u0; u0 = u; u = temp;
    var temp = v0; v0 = v; v = temp;
    this.diffuse2(u,u0,v,v0, dt);
    this.project(u, v, u0, v0);
    var temp = u0; u0 = u; u = temp; 
    var temp = v0; v0 = v; v = temp;
    this.advect(1, u, u0, u0, v0, dt);
    this.advect(2, v, v0, u0, v0, dt);
    this.project(u, v, u0, v0);
    for (var i = 0; i < this.size; i++) {
    u[i] *= this.damp;
    v[i] *= this.damp;
  }
}

FluidField.prototype.setVelocity = function(x, y, xv, yv) {
  this.u[(x + 1) + (y + 1) * this.rowSize] = xv;
  this.v[(x + 1) + (y + 1) * this.rowSize] = yv;
}

FluidField.prototype.getXVelocity = function(x, y) {
  if (this.interpolate) {
    var x_ = Math.floor(x);
    var y_ = Math.floor(y);
    var center = this.u[(x_ + 1) + (y_ + 1) * this.rowSize]; 
    var right = this.u[(x_ + 2) + (y_ + 1) * this.rowSize];
    var up = this.u[(x_ + 1) + (y_ + 2) * this.rowSize];
    var dx = x - x_;
    var dy = y - y_;
    return (center * dx + (1 - dx) * right + center * dy + (1 - dy) * up) / 2;
  } 
  return this.u[(x + 1) + (y + 1) * this.rowSize];
}

FluidField.prototype.getYVelocity = function(x, y) {
  if (this.interpolate) {
    var x_ = Math.floor(x);
    var y_ = Math.floor(y);
    var center = this.v[(x_ + 1) + (y_ + 1) * this.rowSize]; 
    var right = this.v[(x_ + 2) + (y_ + 1) * this.rowSize];
    var up = this.v[(x_ + 1) + (y_ + 2) * this.rowSize];
    var dx = x - x_;
    var dy = y - y_;
    return (center * dx + (1 - dx) * right + center * dy + (1 - dy) * up) / 2;
  }
  return this.v[(x + 1) + (y + 1) * this.rowSize];
}

FluidField.prototype.update = function() {
  this.velocityStep(this.u, this.v, this.u_prev, this.v_prev, this.dt);
}

function FluidFieldDisplay(field, canvas) {
  this.additive = true;
  this.canvas = canvas;
  this.field = field;
  this.context = this.canvas.getContext("2d");
  this.image = null; // created in renderParticles
}

FluidFieldDisplay.prototype.clear = function() {
  this.context.fillStyle = 'black';
  this.context.fillRect(0, 0, this.canvas.width, this.canvas.height);
}

FluidFieldDisplay.prototype.renderVelocityField = function() {
  this.context.save();
  this.context.lineWidth = 1;
  var wscale = this.canvas.width / this.field.width;
  var hscale = this.canvas.height / this.field.height;
  this.context.strokeStyle = 'white';
  var scale = 20;
  this.context.beginPath();
  for (var x = 0; x < this.field.width; x++) {
    for (var y = 0; y < this.field.height; y++) {
      var vx = this.field.getXVelocity(x, y);
      var vy = this.field.getYVelocity(x, y);
      this.context.moveTo(x * wscale + 0.5 * wscale, y * hscale + 0.5 * hscale);
      this.context.lineTo((x + 0.5 + scale * vx) * wscale,
                          (y + 0.5 + scale * vy) * hscale);        
    }
  }
  this.context.stroke();
  this.context.restore();
}

FluidFieldDisplay.prototype.renderParticles = function(field, px, py, pc, pl) {
  this.image = this.context.createImageData(this.canvas.width, this.canvas.height);
  var wscale = this.canvas.width / this.field.width;
  var hscale = this.canvas.height / this.field.height;
  var n = px.length;
  var h, s, v, r, g, b, j, f, p, q, t;
  for (var i = 0; i < n; i++) {
    var cx = Math.floor(px[i] * wscale);
    var cy = Math.floor(py[i] * hscale);
    var base = 4 * (cy * this.canvas.width + cx);
    h = pc[i];
    s = pl[i] / vz.lifetime;
    v = s * s;
    { // hsv -> rgb, inlined for performance
      j = Math.floor(h * 6);
      f = h * 6 - j;
      p = v * (1 - s);
      q = v * (1 - f * s);
      t = v * (1 - (1 - f) * s);
      switch (j % 6) {
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
      }
    }
    // alpha
    this.image.data[base + 3]  = 255;
    // blending modes
    this.image.data[base    ] += Math.floor(r * 255);
    this.image.data[base + 1] += Math.floor(g * 255);
    this.image.data[base + 2] += Math.floor(b * 255);
    this.image.data[base + 4] += Math.floor(r * 255);
    this.image.data[base + 5] += Math.floor(g * 255);
    this.image.data[base + 6] += Math.floor(b * 255);
    this.image.data[base + 7]  = 255;
  }
  this.context.putImageData(this.image, 0, 0);
}
