   var pt, ind, off,pi,pi2, pi0 = 1200000, su = 33,sv = 9, rnd, irnd,
       leaf = true, age = 9, thr1 = 0, thr2 = 0, fir = 0,
       th1 = 0, fi1 = 0, th2 = 0, fi2 = 0, sclen = .8, scr1 = .6,
       anim = 0, bAnim = false, uAnim, tex, green = true, first = true,
       frames = 0, time;
   var prMatrix, mvMat, mvMatLoc, rotMat, c_w, c_h;

   function webGLStart() {
      initGL();
      var err = "Your browser does not support ";
      var ext = gl.getExtension("OES_element_index_uint");
      if ( !ext ) {alert(err + "OES_element_index_uint extension"); return;}
   
      transl = -21;
      c_w = Math.round(.9*window.innerWidth);  c_h = window.innerHeight - 10;
      canvas.width = c_w;   canvas.height = c_h;
      gl.viewport(0, 0, c_w, c_h);
   
      var prog_show  = gl.createProgram();
      gl.attachShader(prog_show, getShader( gl, "shader-vs-show" ));
      gl.attachShader(prog_show, getShader( gl, "shader-fs-show" ));
      gl.linkProgram(prog_show);
      gl.useProgram(prog_show);
      tex = gl.getUniformLocation(prog_show, "uTS");
      uAnim = gl.getUniformLocation(prog_show, "anim");
   
      var posLocation = gl.getAttribLocation(prog_show, "aPos");
      gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
      gl.vertexAttribPointer(posLocation, 3, gl.FLOAT, false, 32, 0);
      gl.enableVertexAttribArray( posLocation );
      var colLocation = gl.getAttribLocation(prog_show, "aCol");
      gl.vertexAttribPointer(colLocation, 3, gl.FLOAT, false, 32, 12);
      gl.enableVertexAttribArray( colLocation );
      var texLoc = gl.getAttribLocation(prog_show, "aTC");
      gl.vertexAttribPointer(texLoc, 2, gl.FLOAT, false, 32, 24);
      gl.enableVertexAttribArray( texLoc );
   
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, gl.createBuffer());
      prMatrix = new CanvasMatrix4();
      prMatrix.perspective(70, c_w/c_h, .01, 100);
      mvMatrix = new CanvasMatrix4();
      rotMat = new CanvasMatrix4();
      rotMat.makeIdentity();
      rotMat.rotate(-160, 0,1,1); // rotate the tree
      mvMatLoc = gl.getUniformLocation(prog_show,"mvMatrix");
   
      var k = 512,  img = new Uint8Array(4* k*k);
      bark_tex(k, img)
      var texture = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, texture);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, k, k, 0, gl.RGBA, gl.UNSIGNED_BYTE, img)
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
      gl.generateMipmap(gl.TEXTURE_2D);
   
      var k = 2*64,  img = new Uint8Array(4* k*k);
      leaf_tex3(k, img)
      gl.activeTexture(gl.TEXTURE1);
      var texture = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, texture);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, k, k, 0, gl.RGBA, gl.UNSIGNED_BYTE, img)
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
      gl.generateMipmap(gl.TEXTURE_2D);
   
      gl.enable(gl.DEPTH_TEST);
      gl.clearColor(.5, .7, 1., 1);
   
      rnd = new Float32Array(10000)
      rand()
      pt = new Float32Array(2700000);  ind = new Uint32Array(1400000);
      sin = new Float32Array(2*su-1);  cos = new Float32Array(2*su-1);
      for(var i = 0; i < 2*su-1; i++ ){
        var fi = Math.PI*i/(su-1);
        sin[i] = Math.sin(fi);  cos[i] = Math.cos(fi);
      }
      tree()
      time = new Date().getTime();
      setInterval(fr, 500);
      animate()
   
     canvas.resize = function (){
       c_w = Math.round(.9*window.innerWidth);  c_h = window.innerHeight - 10;
       canvas.width = c_w;   canvas.height = c_h;
       gl.viewport(0, 0, c_w, c_h);
       prMatrix.makeIdentity();
       prMatrix.perspective(70, c_w/c_h, .01, 100);
       drawScene();
     }
   }
   function animate(){
     if(!bAnim) return
     drawScene();
     anim += .01
     requestAnimationFrame(animate);
   }
   function fr(){
     var ti = new Date().getTime()
     document.getElementById("framerate").value =
       Math.round(1000*frames/(ti - time))
     frames = 0;  time = ti
   }
   function tree(){
      first = true
      off = 0; pi = 0; pi2 = pi0;  irnd = 0
      var b1 = new Float32Array([0,1,0, 0,0,1, 1,0,0, 0,0,0, 0,0]);
      var rg = .1+.05*age
      // console.log('age', age)
      base(b1, 0, 1.5*rg,1.2*rg,.3*rg, .2,1, 1) // to remove the trunk base 
      twig(b1, .1,0, 1.2*rg,rg, 1,1, 2,1)
      branch(age, b1, rg,1.5, 4,[0,1,1,0,0,0,0,0,0], age)
   //alert(pi)
      var t = off/8, i = 0,  r = 1.5; // ground
      var gr = [-r,-r, r,-r, -r,r, r,r];
      for(var j = 0; j < 4; j++ ){
        pt[off++] = gr[i++];  pt[off++] = 0;  pt[off++] = gr[i++];
        pt[off++] = 0; pt[off++] = .6;  pt[off++] = 0;
        pt[off++] = 0; pt[off++] = 0; }
      ind[pi2++] = t++;  ind[pi2++] = t;  ind[pi2++] = t+1;
      ind[pi2++] = t++; ind[pi2++] = t++; ind[pi2++] = t;
   
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, ind, gl.STATIC_DRAW);
      gl.bufferData(gl.ARRAY_BUFFER, pt, gl.STATIC_DRAW);
      drawScene()
   }
   function branch(it, b, sr,sb, stu,div, age){
    // console.log('it', it)
    var b2 = new Float32Array(b),  sr2 = sr*3, age2 = age;
    // console.log('age2', age2)
      rami(b,b2, -thr1,thr2,fir, sr,.8*sr,.8*sr, sr2,sr2,sr2, .3, stu)
      var di = div[it], sr1 = scr1*sr, sr2 = .2*sr
      if(first) twig(b, th1*.7,fi1, .8*sr,sr1, sb,sb, stu,di)
      else twig(b, th1,fi1, .8*sr,sr1, sb,sb, stu,di)
      twig(b2, th2,fi2, .8*sr,sr2, sb,sb, stu,di)
      if(di) stu *= 2
      it--
      if(it == 0){
        if(leaf){ 
          leaves3(b, sr1, green, age2); 
          // leaves3(b2, sr2, green);
        }
        return }
      if(first){  first = false
        branch(it, b, sr1, sb*sclen*.6, stu,div, age2)
      }
      else branch(it, b, sr1, sb*sclen, stu,div, age2)
      // branch(it, b2, sr2, sb*sclen, stu,div) // to remove the second branch
   }
   function rand(){
     for(var i = 0; i < 10000; i++ ) rnd[i] = Math.random()
   }
   function drawScene(){
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      rotMat.rotate(xRot/5, 1,0,0);  rotMat.rotate(yRot/5, 0,1,0);
      rotMat.rotate(zRot, 0,0,1);
      yRot = xRot = zRot = 0;
      mvMatrix.makeIdentity();
      mvMatrix.translate(0, -12, 0);
      mvMatrix.multRight( rotMat );
      mvMatrix.translate(-1, 0, transl);
      mvMatrix.multRight( prMatrix );
      gl.uniformMatrix4fv( mvMatLoc, false, new Float32Array(mvMatrix.getAsArray()) );
      gl.uniform1f(uAnim, .1*Math.sin(anim) + .05*Math.sin(2.7*anim));
      gl.uniform1i(tex, 0);
      gl.enable(gl.CULL_FACE);
      gl.drawElements(gl.TRIANGLES, pi, gl.UNSIGNED_INT, 0);
      gl.disable(gl.CULL_FACE);
      gl.uniform1i(tex, 1);
      gl.drawElements(gl.TRIANGLES, pi2 - pi0, gl.UNSIGNED_INT, 4*pi0);
      frames++
   }
   var guiFork = function() { // dat.GUI
     this.animation = bAnim
     this.green = green
     this.age = age
     this.leaves = leaf
     this.bend1 = thr1
     this.bend2 = thr2
     this.rot = fir
     this.length = sclen
   };
   var guiTwig1 = function() {
     this.bend = th1
     this.rot = fi1
     this.scale_r1 = scr1
   };
   var guiTwig2 = function() {
     this.bend = th2
     this.rot = fi2
   };
   function start() {
     webGLStart()
     var gui = new dat.GUI();
     var tFork = new guiFork();
     gui.add(tFork, 'animation').onChange(function(v){bAnim=v;animate()})
     gui.add(tFork, 'age', 1, 9).step(1).onChange(function(v){age=v;tree()})
     gui.add(tFork, 'green').onChange(function(v){green=v;tree()})
     gui.add(tFork, 'leaves').onChange(function(v){leaf=v;tree()})
     gui.add(tFork, 'length', .5, 1).onChange(function(v){sclen=v;tree()})
     var fFork = gui.addFolder('Fork');
     fFork.add(tFork, 'bend1', .1, 1.5).onChange(function(v){thr1=v;tree()})
     fFork.add(tFork, 'bend2', .1, 1.5).onChange(function(v){thr2=v;tree()})
     fFork.add(tFork, 'rot', -9,9).step(1).onChange(function(v){fir=v;tree()})
   //  fFork.open();
     var tTwig1 = new guiTwig1();
     var fTwig1 = gui.addFolder('Twig 1');
     fTwig1.add(tTwig1, 'bend', -1.5, 1.5).onChange(function(v){th1=v;tree()})
     fTwig1.add(tTwig1, 'rot', -9,9).step(1).onChange(function(v){fi1=v;tree()})
     fTwig1.add(tTwig1, 'scale_r1', .4, .7).onChange(function(v){scr1=v;tree()})
     var tTwig2 = new guiTwig2();
     var fTwig2 = gui.addFolder('Twig 2');
     fTwig2.add(tTwig2, 'bend', -1.5, 1.5).onChange(function(v){th2=v;tree()})
     fTwig2.add(tTwig2, 'rot', -9,9).step(1).onChange(function(v){fi2=v;tree()})
   //  gui.close();
   }
