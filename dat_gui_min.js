/**
 * dat-gui JavaScript Controller Library
 * http://code.google.com/p/dat-gui
 *
 * Copyright 2011 Data Arts Team, Google Creative Lab
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 */
var dat = dat || {};
dat.gui = dat.gui || {},
dat.utils = dat.utils || {},
dat.controllers = dat.controllers || {},
dat.dom = dat.dom || {},
dat.color = dat.color || {},
dat.utils.css = function() {
    return {
        load: function(e, t) {
            t = t || document;
            var n = t.createElement("link");
            n.type = "text/css",
            n.rel = "stylesheet",
            n.href = e,
            t.getElementsByTagName("head")[0].appendChild(n)
        },
        inject: function(e, t) {
            t = t || document;
            var n = document.createElement("style");
            n.type = "text/css",
            n.innerHTML = e,
            t.getElementsByTagName("head")[0].appendChild(n)
        }
    }
}(),
dat.utils.common = function() {
    var e = Array.prototype.forEach
      , t = Array.prototype.slice;
    return {
        BREAK: {},
        extend: function(e) {
            return this.each(t.call(arguments, 1), function(t) {
                for (var n in t)
                    this.isUndefined(t[n]) || (e[n] = t[n])
            }, this),
            e
        },
        defaults: function(e) {
            return this.each(t.call(arguments, 1), function(t) {
                for (var n in t)
                    this.isUndefined(e[n]) && (e[n] = t[n])
            }, this),
            e
        },
        compose: function() {
            var e = t.call(arguments);
            return function() {
                for (var n = t.call(arguments), o = e.length - 1; o >= 0; o--)
                    n = [e[o].apply(this, n)];
                return n[0]
            }
        },
        each: function(t, n, o) {
            if (t)
                if (e && t.forEach && t.forEach === e)
                    t.forEach(n, o);
                else if (t.length === t.length + 0) {
                    for (var i = 0, r = t.length; r > i; i++)
                        if (i in t && n.call(o, t[i], i) === this.BREAK)
                            return
                } else
                    for (var i in t)
                        if (n.call(o, t[i], i) === this.BREAK)
                            return
        },
        defer: function(e) {
            setTimeout(e, 0)
        },
        toArray: function(e) {
            return e.toArray ? e.toArray() : t.call(e)
        },
        isUndefined: function(e) {
            return void 0 === e
        },
        isNull: function(e) {
            return null === e
        },
        isNaN: function(e) {
            return e !== e
        },
        isArray: Array.isArray || function(e) {
            return e.constructor === Array
        }
        ,
        isObject: function(e) {
            return e === Object(e)
        },
        isNumber: function(e) {
            return e === e + 0
        },
        isString: function(e) {
            return e === e + ""
        },
        isBoolean: function(e) {
            return e === !1 || e === !0
        },
        isFunction: function(e) {
            return "[object Function]" === Object.prototype.toString.call(e)
        }
    }
}(),
dat.controllers.Controller = function(e) {
    var t = function(e, t) {
        this.initialValue = e[t],
        this.domElement = document.createElement("div"),
        this.object = e,
        this.property = t,
        this.__onChange = void 0,
        this.__onFinishChange = void 0
    };
    return e.extend(t.prototype, {
        onChange: function(e) {
            return this.__onChange = e,
            this
        },
        onFinishChange: function(e) {
            return this.__onFinishChange = e,
            this
        },
        setValue: function(e) {
            return this.object[this.property] = e,
            this.__onChange && this.__onChange.call(this, e),
            this.updateDisplay(),
            this
        },
        getValue: function() {
            return this.object[this.property]
        },
        updateDisplay: function() {
            return this
        },
        isModified: function() {
            return this.initialValue !== this.getValue()
        }
    }),
    t
}(dat.utils.common),
dat.dom.dom = function(e) {
    function t(t) {
        if ("0" === t || e.isUndefined(t))
            return 0;
        var n = t.match(i);
        return e.isNull(n) ? 0 : parseFloat(n[1])
    }
    var n = {
        HTMLEvents: ["change"],
        MouseEvents: ["click", "mousemove", "mousedown", "mouseup", "mouseover"],
        KeyboardEvents: ["keydown"]
    }
      , o = {};
    e.each(n, function(t, n) {
        e.each(t, function(e) {
            o[e] = n
        })
    });
    var i = /(\d+(\.\d+)?)px/
      , r = {
        makeSelectable: function(e, t) {
            void 0 !== e && void 0 !== e.style && (e.onselectstart = t ? function() {
                return !1
            }
            : function() {}
            ,
            e.style.MozUserSelect = t ? "auto" : "none",
            e.style.KhtmlUserSelect = t ? "auto" : "none",
            e.unselectable = t ? "on" : "off")
        },
        makeFullscreen: function(t, n, o) {
            e.isUndefined(n) && (n = !0),
            e.isUndefined(o) && (o = !0),
            t.style.position = "absolute",
            n && (t.style.left = 0,
            t.style.right = 0),
            o && (t.style.top = 0,
            t.style.bottom = 0)
        },
        fakeEvent: function(t, n, i, r) {
            i = i || {};
            var s = o[n];
            if (!s)
                throw new Error("Event type " + n + " not supported.");
            var a = document.createEvent(s);
            switch (s) {
            case "MouseEvents":
                var l = i.x || i.clientX || 0
                  , d = i.y || i.clientY || 0;
                a.initMouseEvent(n, i.bubbles || !1, i.cancelable || !0, window, i.clickCount || 1, 0, 0, l, d, !1, !1, !1, !1, 0, null);
                break;
            case "KeyboardEvents":
                var c = a.initKeyboardEvent || a.initKeyEvent;
                e.defaults(i, {
                    cancelable: !0,
                    ctrlKey: !1,
                    altKey: !1,
                    shiftKey: !1,
                    metaKey: !1,
                    keyCode: void 0,
                    charCode: void 0
                }),
                c(n, i.bubbles || !1, i.cancelable, window, i.ctrlKey, i.altKey, i.shiftKey, i.metaKey, i.keyCode, i.charCode);
                break;
            default:
                a.initEvent(n, i.bubbles || !1, i.cancelable || !0)
            }
            e.defaults(a, r),
            t.dispatchEvent(a)
        },
        bind: function(e, t, n, o) {
            return o = o || !1,
            e.addEventListener ? e.addEventListener(t, n, o) : e.attachEvent && e.attachEvent("on" + t, n),
            r
        },
        unbind: function(e, t, n, o) {
            return o = o || !1,
            e.removeEventListener ? e.removeEventListener(t, n, o) : e.detachEvent && e.detachEvent("on" + t, n),
            r
        },
        addClass: function(e, t) {
            if (void 0 === e.className)
                e.className = t;
            else if (e.className !== t) {
                var n = e.className.split(/ +/);
                -1 == n.indexOf(t) && (n.push(t),
                e.className = n.join(" ").replace(/^\s+/, "").replace(/\s+$/, ""))
            }
            return r
        },
        removeClass: function(e, t) {
            if (t)
                if (void 0 === e.className)
                    ;
                else if (e.className === t)
                    e.removeAttribute("class");
                else {
                    var n = e.className.split(/ +/)
                      , o = n.indexOf(t);
                    -1 != o && (n.splice(o, 1),
                    e.className = n.join(" "))
                }
            else
                e.className = void 0;
            return r
        },
        hasClass: function(e, t) {
            return new RegExp("(?:^|\\s+)" + t + "(?:\\s+|$)").test(e.className) || !1
        },
        getWidth: function(e) {
            var n = getComputedStyle(e);
            return t(n["border-left-width"]) + t(n["border-right-width"]) + t(n["padding-left"]) + t(n["padding-right"]) + t(n.width)
        },
        getHeight: function(e) {
            var n = getComputedStyle(e);
            return t(n["border-top-width"]) + t(n["border-bottom-width"]) + t(n["padding-top"]) + t(n["padding-bottom"]) + t(n.height)
        },
        getOffset: function(e) {
            var t = {
                left: 0,
                top: 0
            };
            if (e.offsetParent)
                do
                    t.left += e.offsetLeft,
                    t.top += e.offsetTop;
                while (e = e.offsetParent);
            return t
        },
        isActive: function(e) {
            return e === document.activeElement && (e.type || e.href)
        }
    };
    return r
}(dat.utils.common),
dat.controllers.OptionController = function(e, t, n) {
    var o = function(e, i, r) {
        o.superclass.call(this, e, i);
        var s = this;
        if (this.__select = document.createElement("select"),
        n.isArray(r)) {
            var a = {};
            n.each(r, function(e) {
                a[e] = e
            }),
            r = a
        }
        n.each(r, function(e, t) {
            var n = document.createElement("option");
            n.innerHTML = t,
            n.setAttribute("value", e),
            s.__select.appendChild(n)
        }),
        this.updateDisplay(),
        t.bind(this.__select, "change", function() {
            var e = this.options[this.selectedIndex].value;
            s.setValue(e)
        }),
        this.domElement.appendChild(this.__select)
    };
    return o.superclass = e,
    n.extend(o.prototype, e.prototype, {
        setValue: function(e) {
            var t = o.superclass.prototype.setValue.call(this, e);
            return this.__onFinishChange && this.__onFinishChange.call(this, this.getValue()),
            t
        },
        updateDisplay: function() {
            return this.__select.value = this.getValue(),
            o.superclass.prototype.updateDisplay.call(this)
        }
    }),
    o
}(dat.controllers.Controller, dat.dom.dom, dat.utils.common),
dat.controllers.NumberController = function(e, t) {
    function n(e) {
        return e = e.toString(),
        e.indexOf(".") > -1 ? e.length - e.indexOf(".") - 1 : 0
    }
    var o = function(e, i, r) {
        o.superclass.call(this, e, i),
        r = r || {},
        this.__min = r.min,
        this.__max = r.max,
        this.__step = r.step,
        this.__impliedStep = t.isUndefined(this.__step) ? 0 == this.initialValue ? 1 : Math.pow(10, Math.floor(Math.log(Math.abs(this.initialValue)) / Math.LN10)) / 10 : this.__step,
        this.__precision = n(this.__impliedStep)
    };
    return o.superclass = e,
    t.extend(o.prototype, e.prototype, {
        setValue: function(e) {
            return void 0 !== this.__min && e < this.__min ? e = this.__min : void 0 !== this.__max && e > this.__max && (e = this.__max),
            void 0 !== this.__step && e % this.__step != 0 && (e = Math.round(e / this.__step) * this.__step),
            o.superclass.prototype.setValue.call(this, e)
        },
        min: function(e) {
            return this.__min = e,
            this
        },
        max: function(e) {
            return this.__max = e,
            this
        },
        step: function(e) {
            return this.__step = e,
            this.__impliedStep = e,
            this.__precision = n(e),
            this
        }
    }),
    o
}(dat.controllers.Controller, dat.utils.common),
dat.controllers.NumberControllerBox = function(e, t, n) {
    function o(e, t) {
        var n = Math.pow(10, t);
        return Math.round(e * n) / n
    }
    var i = function(e, o, r) {
        function s() {
            var e = parseFloat(h.__input.value);
            n.isNaN(e) || h.setValue(e)
        }
        function a() {
            s(),
            h.__onFinishChange && h.__onFinishChange.call(h, h.getValue())
        }
        function l(e) {
            t.bind(window, "mousemove", d),
            t.bind(window, "mouseup", c),
            u = e.clientY
        }
        function d(e) {
            var t = u - e.clientY;
            h.setValue(h.getValue() + t * h.__impliedStep),
            u = e.clientY
        }
        function c() {
            t.unbind(window, "mousemove", d),
            t.unbind(window, "mouseup", c)
        }
        this.__truncationSuspended = !1,
        i.superclass.call(this, e, o, r);
        var u, h = this;
        this.__input = document.createElement("input"),
        this.__input.setAttribute("type", "text"),
        t.bind(this.__input, "change", s),
        t.bind(this.__input, "blur", a),
        t.bind(this.__input, "mousedown", l),
        t.bind(this.__input, "keydown", function(e) {
            13 === e.keyCode && (h.__truncationSuspended = !0,
            this.blur(),
            h.__truncationSuspended = !1)
        }),
        this.updateDisplay(),
        this.domElement.appendChild(this.__input)
    };
    return i.superclass = e,
    n.extend(i.prototype, e.prototype, {
        updateDisplay: function() {
            return this.__input.value = this.__truncationSuspended ? this.getValue() : o(this.getValue(), this.__precision),
            i.superclass.prototype.updateDisplay.call(this)
        }
    }),
    i
}(dat.controllers.NumberController, dat.dom.dom, dat.utils.common),
dat.controllers.NumberControllerSlider = function(e, t, n, o, i) {
    function r(e, t, n, o, i) {
        return o + (i - o) * ((e - t) / (n - t))
    }
    var s = function(e, n, o, i, a) {
        function l(e) {
            t.bind(window, "mousemove", d),
            t.bind(window, "mouseup", c),
            d(e)
        }
        function d(e) {
            e.preventDefault();
            var n = t.getOffset(u.__background)
              , o = t.getWidth(u.__background);
            return u.setValue(r(e.clientX, n.left, n.left + o, u.__min, u.__max)),
            !1
        }
        function c() {
            t.unbind(window, "mousemove", d),
            t.unbind(window, "mouseup", c),
            u.__onFinishChange && u.__onFinishChange.call(u, u.getValue())
        }
        s.superclass.call(this, e, n, {
            min: o,
            max: i,
            step: a
        });
        var u = this;
        this.__background = document.createElement("div"),
        this.__foreground = document.createElement("div"),
        t.bind(this.__background, "mousedown", l),
        t.addClass(this.__background, "slider"),
        t.addClass(this.__foreground, "slider-fg"),
        this.updateDisplay(),
        this.__background.appendChild(this.__foreground),
        this.domElement.appendChild(this.__background)
    };
    return s.superclass = e,
    s.useDefaultStyles = function() {
        n.inject(i)
    }
    ,
    o.extend(s.prototype, e.prototype, {
        updateDisplay: function() {
            var e = (this.getValue() - this.__min) / (this.__max - this.__min);
            return this.__foreground.style.width = 100 * e + "%",
            s.superclass.prototype.updateDisplay.call(this)
        }
    }),
    s
}(dat.controllers.NumberController, dat.dom.dom, dat.utils.css, dat.utils.common, "/**\n * dat-gui JavaScript Controller Library\n * http://code.google.com/p/dat-gui\n *\n * Copyright 2011 Data Arts Team, Google Creative Lab\n *\n * Licensed under the Apache License, Version 2.0 (the \"License\");\n * you may not use this file except in compliance with the License.\n * You may obtain a copy of the License at\n *\n * http://www.apache.org/licenses/LICENSE-2.0\n */\n\n.slider {\n  box-shadow: inset 0 2px 4px rgba(0,0,0,0.15);\n  height: 1em;\n  border-radius: 1em;\n  background-color: #eee;\n  padding: 0 0.5em;\n  overflow: hidden;\n}\n\n.slider-fg {\n  padding: 1px 0 2px 0;\n  background-color: #aaa;\n  height: 1em;\n  margin-left: -0.5em;\n  padding-right: 0.5em;\n  border-radius: 1em 0 0 1em;\n}\n\n.slider-fg:after {\n  display: inline-block;\n  border-radius: 1em;\n  background-color: #fff;\n  border:  1px solid #aaa;\n  content: '';\n  float: right;\n  margin-right: -1em;\n  margin-top: -1px;\n  height: 0.9em;\n  width: 0.9em;\n}"),
dat.controllers.FunctionController = function(e, t, n) {
    var o = function(e, n, i) {
        o.superclass.call(this, e, n);
        var r = this;
        this.__button = document.createElement("div"),
        this.__button.innerHTML = void 0 === i ? "Fire" : i,
        t.bind(this.__button, "click", function(e) {
            return e.preventDefault(),
            r.fire(),
            !1
        }),
        t.addClass(this.__button, "button"),
        this.domElement.appendChild(this.__button)
    };
    return o.superclass = e,
    n.extend(o.prototype, e.prototype, {
        fire: function() {
            this.__onChange && this.__onChange.call(this),
            this.getValue().call(this.object),
            this.__onFinishChange && this.__onFinishChange.call(this, this.getValue())
        }
    }),
    o
}(dat.controllers.Controller, dat.dom.dom, dat.utils.common),
dat.controllers.BooleanController = function(e, t, n) {
    var o = function(e, n) {
        function i() {
            r.setValue(!r.__prev)
        }
        o.superclass.call(this, e, n);
        var r = this;
        this.__prev = this.getValue(),
        this.__checkbox = document.createElement("input"),
        this.__checkbox.setAttribute("type", "checkbox"),
        t.bind(this.__checkbox, "change", i, !1),
        this.domElement.appendChild(this.__checkbox),
        this.updateDisplay()
    };
    return o.superclass = e,
    n.extend(o.prototype, e.prototype, {
        setValue: function(e) {
            var t = o.superclass.prototype.setValue.call(this, e);
            return this.__onFinishChange && this.__onFinishChange.call(this, this.getValue()),
            this.__prev = this.getValue(),
            t
        },
        updateDisplay: function() {
            return this.getValue() === !0 ? (this.__checkbox.setAttribute("checked", "checked"),
            this.__checkbox.checked = !0) : this.__checkbox.checked = !1,
            o.superclass.prototype.updateDisplay.call(this)
        }
    }),
    o
}(dat.controllers.Controller, dat.dom.dom, dat.utils.common),
dat.color.toString = function(e) {
    return function(t) {
        if (1 == t.a || e.isUndefined(t.a)) {
            for (var n = t.hex.toString(16); n.length < 6; )
                n = "0" + n;
            return "#" + n
        }
        return "rgba(" + Math.round(t.r) + "," + Math.round(t.g) + "," + Math.round(t.b) + "," + t.a + ")"
    }
}(dat.utils.common),
dat.color.interpret = function(e, t) {
    var n, o, i = function() {
        o = !1;
        var e = arguments.length > 1 ? t.toArray(arguments) : arguments[0];
        return t.each(r, function(i) {
            return i.litmus(e) ? (t.each(i.conversions, function(i, r) {
                return n = i.read(e),
                o === !1 && n !== !1 ? (o = n,
                n.conversionName = r,
                n.conversion = i,
                t.BREAK) : void 0
            }),
            t.BREAK) : void 0
        }),
        o
    }, r = [{
        litmus: t.isString,
        conversions: {
            THREE_CHAR_HEX: {
                read: function(e) {
                    var t = e.match(/^#([A-F0-9])([A-F0-9])([A-F0-9])$/i);
                    return null === t ? !1 : {
                        space: "HEX",
                        hex: parseInt("0x" + t[1].toString() + t[1].toString() + t[2].toString() + t[2].toString() + t[3].toString() + t[3].toString())
                    }
                },
                write: e
            },
            SIX_CHAR_HEX: {
                read: function(e) {
                    var t = e.match(/^#([A-F0-9]{6})$/i);
                    return null === t ? !1 : {
                        space: "HEX",
                        hex: parseInt("0x" + t[1].toString())
                    }
                },
                write: e
            },
            CSS_RGB: {
                read: function(e) {
                    var t = e.match(/^rgb\(\s*(.+)\s*,\s*(.+)\s*,\s*(.+)\s*\)/);
                    return null === t ? !1 : {
                        space: "RGB",
                        r: parseFloat(t[1]),
                        g: parseFloat(t[2]),
                        b: parseFloat(t[3])
                    }
                },
                write: e
            },
            CSS_RGBA: {
                read: function(e) {
                    var t = e.match(/^rgba\(\s*(.+)\s*,\s*(.+)\s*,\s*(.+)\s*\,\s*(.+)\s*\)/);
                    return null === t ? !1 : {
                        space: "RGB",
                        r: parseFloat(t[1]),
                        g: parseFloat(t[2]),
                        b: parseFloat(t[3]),
                        a: parseFloat(t[4])
                    }
                },
                write: e
            }
        }
    }, {
        litmus: t.isNumber,
        conversions: {
            HEX: {
                read: function(e) {
                    return {
                        space: "HEX",
                        hex: e,
                        conversionName: "HEX"
                    }
                },
                write: function(e) {
                    return e.hex
                }
            }
        }
    }, {
        litmus: t.isArray,
        conversions: {
            RGB_ARRAY: {
                read: function(e) {
                    return 3 != e.length ? !1 : {
                        space: "RGB",
                        r: e[0],
                        g: e[1],
                        b: e[2]
                    }
                },
                write: function(e) {
                    return [e.r, e.g, e.b]
                }
            },
            RGBA_ARRAY: {
                read: function(e) {
                    return 4 != e.length ? !1 : {
                        space: "RGB",
                        r: e[0],
                        g: e[1],
                        b: e[2],
                        a: e[3]
                    }
                },
                write: function(e) {
                    return [e.r, e.g, e.b, e.a]
                }
            }
        }
    }, {
        litmus: t.isObject,
        conversions: {
            RGBA_OBJ: {
                read: function(e) {
                    return t.isNumber(e.r) && t.isNumber(e.g) && t.isNumber(e.b) && t.isNumber(e.a) ? {
                        space: "RGB",
                        r: e.r,
                        g: e.g,
                        b: e.b,
                        a: e.a
                    } : !1
                },
                write: function(e) {
                    return {
                        r: e.r,
                        g: e.g,
                        b: e.b,
                        a: e.a
                    }
                }
            },
            RGB_OBJ: {
                read: function(e) {
                    return t.isNumber(e.r) && t.isNumber(e.g) && t.isNumber(e.b) ? {
                        space: "RGB",
                        r: e.r,
                        g: e.g,
                        b: e.b
                    } : !1
                },
                write: function(e) {
                    return {
                        r: e.r,
                        g: e.g,
                        b: e.b
                    }
                }
            },
            HSVA_OBJ: {
                read: function(e) {
                    return t.isNumber(e.h) && t.isNumber(e.s) && t.isNumber(e.v) && t.isNumber(e.a) ? {
                        space: "HSV",
                        h: e.h,
                        s: e.s,
                        v: e.v,
                        a: e.a
                    } : !1
                },
                write: function(e) {
                    return {
                        h: e.h,
                        s: e.s,
                        v: e.v,
                        a: e.a
                    }
                }
            },
            HSV_OBJ: {
                read: function(e) {
                    return t.isNumber(e.h) && t.isNumber(e.s) && t.isNumber(e.v) ? {
                        space: "HSV",
                        h: e.h,
                        s: e.s,
                        v: e.v
                    } : !1
                },
                write: function(e) {
                    return {
                        h: e.h,
                        s: e.s,
                        v: e.v
                    }
                }
            }
        }
    }];
    return i
}(dat.color.toString, dat.utils.common),
dat.GUI = dat.gui.GUI = function(e, t, n, o, i, r, s, a, l, d, c, u, h, p, _) {
    function f(e, t, n, r) {
        if (void 0 === t[n])
            throw new Error("Object " + t + ' has no property "' + n + '"');
        var s;
        if (r.color)
            s = new c(t,n);
        else {
            var a = [t, n].concat(r.factoryArgs);
            s = o.apply(e, a)
        }
        r.before instanceof i && (r.before = r.before.__li),
        b(e, s),
        p.addClass(s.domElement, "c");
        var l = document.createElement("span");
        p.addClass(l, "property-name"),
        l.innerHTML = s.property;
        var d = document.createElement("div");
        d.appendChild(l),
        d.appendChild(s.domElement);
        var u = m(e, d, r.before);
        return p.addClass(u, H.CLASS_CONTROLLER_ROW),
        p.addClass(u, typeof s.getValue()),
        g(e, u, s),
        e.__controllers.push(s),
        s
    }
    function m(e, t, n) {
        var o = document.createElement("li");
        return t && o.appendChild(t),
        n ? e.__ul.insertBefore(o, params.before) : e.__ul.appendChild(o),
        e.onResize(),
        o
    }
    function g(e, t, n) {
        if (n.__li = t,
        n.__gui = e,
        _.extend(n, {
            options: function(t) {
                return arguments.length > 1 ? (n.remove(),
                f(e, n.object, n.property, {
                    before: n.__li.nextElementSibling,
                    factoryArgs: [_.toArray(arguments)]
                })) : _.isArray(t) || _.isObject(t) ? (n.remove(),
                f(e, n.object, n.property, {
                    before: n.__li.nextElementSibling,
                    factoryArgs: [t]
                })) : void 0
            },
            name: function(e) {
                return n.__li.firstElementChild.firstElementChild.innerHTML = e,
                n
            },
            listen: function() {
                return n.__gui.listen(n),
                n
            },
            remove: function() {
                return n.__gui.remove(n),
                n
            }
        }),
        n instanceof l) {
            var o = new a(n.object,n.property,{
                min: n.__min,
                max: n.__max,
                step: n.__step
            });
            _.each(["updateDisplay", "onChange", "onFinishChange"], function(e) {
                var t = n[e]
                  , i = o[e];
                n[e] = o[e] = function() {
                    var e = Array.prototype.slice.call(arguments);
                    return t.apply(n, e),
                    i.apply(o, e)
                }
            }),
            p.addClass(t, "has-slider"),
            n.domElement.insertBefore(o.domElement, n.domElement.firstElementChild)
        } else if (n instanceof a) {
            var i = function(t) {
                return _.isNumber(n.__min) && _.isNumber(n.__max) ? (n.remove(),
                f(e, n.object, n.property, {
                    before: n.__li.nextElementSibling,
                    factoryArgs: [n.__min, n.__max, n.__step]
                })) : t
            };
            n.min = _.compose(i, n.min),
            n.max = _.compose(i, n.max)
        } else
            n instanceof r ? (p.bind(t, "click", function() {
                p.fakeEvent(n.__checkbox, "click")
            }),
            p.bind(n.__checkbox, "click", function(e) {
                e.stopPropagation()
            })) : n instanceof s ? (p.bind(t, "click", function() {
                p.fakeEvent(n.__button, "click")
            }),
            p.bind(t, "mouseover", function() {
                p.addClass(n.__button, "hover")
            }),
            p.bind(t, "mouseout", function() {
                p.removeClass(n.__button, "hover")
            })) : n instanceof c && (p.addClass(t, "color"),
            n.updateDisplay = _.compose(function(e) {
                return t.style.borderLeftColor = n.__color.toString(),
                e
            }, n.updateDisplay),
            n.updateDisplay());
        n.setValue = _.compose(function(t) {
            return e.getRoot().__preset_select && n.isModified() && k(e.getRoot(), !0),
            t
        }, n.setValue)
    }
    function b(e, t) {
        var n = e.getRoot()
          , o = n.__rememberedObjects.indexOf(t.object);
        if (-1 != o) {
            var i = n.__rememberedObjectIndecesToControllers[o];
            if (void 0 === i && (i = {},
            n.__rememberedObjectIndecesToControllers[o] = i),
            i[t.property] = t,
            n.load && n.load.remembered) {
                var r, s = n.load.remembered;
                if (s[e.preset])
                    r = s[e.preset];
                else {
                    if (!s[B])
                        return;
                    r = s[B]
                }
                if (r[o] && void 0 !== r[o][t.property]) {
                    var a = r[o][t.property];
                    t.initialValue = a,
                    t.setValue(a)
                }
            }
        }
    }
    function v(e, t) {
        return document.location.href + "." + t
    }
    function y(e) {
        function t() {
            d.style.display = e.useLocalStorage ? "block" : "none"
        }
        var n = e.__save_row = document.createElement("li");
        p.addClass(e.domElement, "has-save"),
        e.__ul.insertBefore(n, e.__ul.firstChild),
        p.addClass(n, "save-row");
        var o = document.createElement("span");
        o.innerHTML = "&nbsp;",
        p.addClass(o, "button gears");
        var i = document.createElement("span");
        i.innerHTML = "Save",
        p.addClass(i, "button"),
        p.addClass(i, "save");
        var r = document.createElement("span");
        r.innerHTML = "New",
        p.addClass(r, "button"),
        p.addClass(r, "save-as");
        var s = document.createElement("span");
        s.innerHTML = "Revert",
        p.addClass(s, "button"),
        p.addClass(s, "revert");
        var a = e.__preset_select = document.createElement("select");
        if (e.load && e.load.remembered ? _.each(e.load.remembered, function(t, n) {
            E(e, n, n == e.preset)
        }) : E(e, B, !1),
        p.bind(a, "change", function() {
            for (var t = 0; t < e.__preset_select.length; t++)
                e.__preset_select[t].innerHTML = e.__preset_select[t].value;
            e.preset = this.value
        }),
        n.appendChild(a),
        n.appendChild(o),
        n.appendChild(i),
        n.appendChild(r),
        n.appendChild(s),
        D) {
            var l = document.getElementById("dg-save-locally")
              , d = document.getElementById("dg-local-explain");
            l.style.display = "block";
            var c = document.getElementById("dg-local-storage");
            "true" === localStorage.getItem(v(e, "isLocal")) && c.setAttribute("checked", "checked"),
            t(),
            p.bind(c, "change", function() {
                e.useLocalStorage = !e.useLocalStorage,
                t()
            })
        }
        var u = document.getElementById("dg-new-constructor");
        p.bind(u, "keydown", function(e) {
            !e.metaKey || 67 !== e.which && 67 != e.keyCode || O.hide()
        }),
        p.bind(o, "click", function() {
            u.innerHTML = JSON.stringify(e.getSaveObject(), void 0, 2),
            O.show(),
            u.focus(),
            u.select()
        }),
        p.bind(i, "click", function() {
            e.save()
        }),
        p.bind(r, "click", function() {
            var t = prompt("Enter a new preset name.");
            t && e.saveAs(t)
        }),
        p.bind(s, "click", function() {
            e.revert()
        })
    }
    function x(e) {
        function t(t) {
            return t.preventDefault(),
            i = t.clientX,
            p.addClass(e.__closeButton, H.CLASS_DRAG),
            p.bind(window, "mousemove", n),
            p.bind(window, "mouseup", o),
            !1
        }
        function n(t) {
            return t.preventDefault(),
            e.width += i - t.clientX,
            e.onResize(),
            i = t.clientX,
            !1
        }
        function o() {
            p.removeClass(e.__closeButton, H.CLASS_DRAG),
            p.unbind(window, "mousemove", n),
            p.unbind(window, "mouseup", o)
        }
        e.__resize_handle = document.createElement("div"),
        _.extend(e.__resize_handle.style, {
            width: "6px",
            marginLeft: "-3px",
            height: "200px",
            cursor: "ew-resize",
            position: "absolute"
        });
        var i;
        p.bind(e.__resize_handle, "mousedown", t),
        p.bind(e.__closeButton, "mousedown", t),
        e.domElement.insertBefore(e.__resize_handle, e.domElement.firstElementChild)
    }
    function w(e, t) {
        e.domElement.style.width = t + "px",
        e.__save_row && e.autoPlace && (e.__save_row.style.width = t + "px"),
        e.__closeButton && (e.__closeButton.style.width = t + "px")
    }
    function C(e, t) {
        var n = {};
        return _.each(e.__rememberedObjects, function(o, i) {
            var r = {}
              , s = e.__rememberedObjectIndecesToControllers[i];
            _.each(s, function(e, n) {
                r[n] = t ? e.initialValue : e.getValue()
            }),
            n[i] = r
        }),
        n
    }
    function E(e, t, n) {
        var o = document.createElement("option");
        o.innerHTML = t,
        o.value = t,
        e.__preset_select.appendChild(o),
        n && (e.__preset_select.selectedIndex = e.__preset_select.length - 1)
    }
    function A(e) {
        for (var t = 0; t < e.__preset_select.length; t++)
            e.__preset_select[t].value == e.preset && (e.__preset_select.selectedIndex = t)
    }
    function k(e, t) {
        var n = e.__preset_select[e.__preset_select.selectedIndex];
        n.innerHTML = t ? n.value + "*" : n.value
    }
    function S(e) {
        0 != e.length && u(function() {
            S(e)
        }),
        _.each(e, function(e) {
            e.updateDisplay()
        })
    }
    e.inject(n);
    var O, T, L = "dg", N = 72, R = 20, B = "Default", D = function() {
        try {
            return "localStorage"in window && null !== window.localStorage
        } catch (e) {
            return !1
        }
    }(), F = !0, V = !1, I = [], H = function(e) {
        function t() {
            var e = n.getRoot();
            e.width += 1,
            _.defer(function() {
                e.width -= 1
            })
        }
        var n = this;
        this.domElement = document.createElement("div"),
        this.__ul = document.createElement("ul"),
        this.domElement.appendChild(this.__ul),
        p.addClass(this.domElement, L),
        this.__folders = {},
        this.__controllers = [],
        this.__rememberedObjects = [],
        this.__rememberedObjectIndecesToControllers = [],
        this.__listening = [],
        e = e || {},
        e = _.defaults(e, {
            autoPlace: !0,
            width: H.DEFAULT_WIDTH
        }),
        e = _.defaults(e, {
            resizable: e.autoPlace,
            hideable: e.autoPlace
        }),
        _.isUndefined(e.load) ? e.load = {
            preset: B
        } : e.preset && (e.load.preset = e.preset),
        _.isUndefined(e.parent) && e.hideable && I.push(this),
        e.resizable = _.isUndefined(e.parent) && e.resizable,
        e.autoPlace && _.isUndefined(e.scrollable) && (e.scrollable = !0);
        var o, i = D && "true" === localStorage.getItem(v(this, "isLocal"));
        if (Object.defineProperties(this, {
            parent: {
                get: function() {
                    return e.parent
                }
            },
            scrollable: {
                get: function() {
                    return e.scrollable
                }
            },
            autoPlace: {
                get: function() {
                    return e.autoPlace
                }
            },
            preset: {
                get: function() {
                    return n.parent ? n.getRoot().preset : e.load.preset
                },
                set: function(t) {
                    n.parent ? n.getRoot().preset = t : e.load.preset = t,
                    A(this),
                    n.revert()
                }
            },
            width: {
                get: function() {
                    return e.width
                },
                set: function(t) {
                    e.width = t,
                    w(n, t)
                }
            },
            name: {
                get: function() {
                    return e.name
                },
                set: function(t) {
                    e.name = t,
                    s && (s.innerHTML = e.name)
                }
            },
            closed: {
                get: function() {
                    return e.closed
                },
                set: function(t) {
                    e.closed = t,
                    e.closed ? p.addClass(n.__ul, H.CLASS_CLOSED) : p.removeClass(n.__ul, H.CLASS_CLOSED),
                    this.onResize(),
                    n.__closeButton && (n.__closeButton.innerHTML = t ? H.TEXT_OPEN : H.TEXT_CLOSED)
                }
            },
            load: {
                get: function() {
                    return e.load
                }
            },
            useLocalStorage: {
                get: function() {
                    return i
                },
                set: function(e) {
                    D && (i = e,
                    e ? p.bind(window, "unload", o) : p.unbind(window, "unload", o),
                    localStorage.setItem(v(n, "isLocal"), e))
                }
            }
        }),
        _.isUndefined(e.parent)) {
            if (e.closed = !1,
            p.addClass(this.domElement, H.CLASS_MAIN),
            p.makeSelectable(this.domElement, !1),
            D && i) {
                n.useLocalStorage = !0;
                var r = localStorage.getItem(v(this, "gui"));
                r && (e.load = JSON.parse(r))
            }
            this.__closeButton = document.createElement("div"),
            this.__closeButton.innerHTML = H.TEXT_CLOSED,
            p.addClass(this.__closeButton, H.CLASS_CLOSE_BUTTON),
            this.domElement.appendChild(this.__closeButton),
            p.bind(this.__closeButton, "click", function() {
                n.closed = !n.closed
            })
        } else {
            void 0 === e.closed && (e.closed = !0);
            var s = document.createTextNode(e.name);
            p.addClass(s, "controller-name");
            var a = m(n, s)
              , l = function(e) {
                return e.preventDefault(),
                n.closed = !n.closed,
                !1
            };
            p.addClass(this.__ul, H.CLASS_CLOSED),
            p.addClass(a, "title"),
            p.bind(a, "click", l),
            e.closed || (this.closed = !1)
        }
        e.autoPlace && (_.isUndefined(e.parent) && (F && (T = document.createElement("div"),
        p.addClass(T, L),
        p.addClass(T, H.CLASS_AUTO_PLACE_CONTAINER),
        document.body.appendChild(T),
        F = !1),
        T.appendChild(this.domElement),
        p.addClass(this.domElement, H.CLASS_AUTO_PLACE)),
        this.parent || w(n, e.width)),
        p.bind(window, "resize", function() {
            n.onResize()
        }),
        p.bind(this.__ul, "webkitTransitionEnd", function() {
            n.onResize()
        }),
        p.bind(this.__ul, "transitionend", function() {
            n.onResize()
        }),
        p.bind(this.__ul, "oTransitionEnd", function() {
            n.onResize()
        }),
        this.onResize(),
        e.resizable && x(this),
        o = function() {
            D && "true" === localStorage.getItem(v(n, "isLocal")) && localStorage.setItem(v(n, "gui"), JSON.stringify(n.getSaveObject()))
        }
        ,
        this.saveToLocalStorageIfPossible = o;
        n.getRoot();
        e.parent || t()
    };
    return H.toggleHide = function() {
        V = !V,
        _.each(I, function(e) {
            e.domElement.style.zIndex = V ? -999 : 999,
            e.domElement.style.opacity = V ? 0 : 1
        })
    }
    ,
    H.CLASS_AUTO_PLACE = "a",
    H.CLASS_AUTO_PLACE_CONTAINER = "ac",
    H.CLASS_MAIN = "main",
    H.CLASS_CONTROLLER_ROW = "cr",
    H.CLASS_TOO_TALL = "taller-than-window",
    H.CLASS_CLOSED = "closed",
    H.CLASS_CLOSE_BUTTON = "close-button",
    H.CLASS_DRAG = "drag",
    H.DEFAULT_WIDTH = 245,
    H.TEXT_CLOSED = "Close Controls",
    H.TEXT_OPEN = "Open Controls",
    p.bind(window, "keydown", function(e) {
        "text" === document.activeElement.type || e.which !== N && e.keyCode != N || H.toggleHide()
    }, !1),
    _.extend(H.prototype, {
        add: function(e, t) {
            return f(this, e, t, {
                factoryArgs: Array.prototype.slice.call(arguments, 2)
            })
        },
        addColor: function(e, t) {
            return f(this, e, t, {
                color: !0
            })
        },
        remove: function(e) {
            this.__ul.removeChild(e.__li),
            this.__controllers.splice(this.__controllers.indexOf(e), 1);
            var t = this;
            _.defer(function() {
                t.onResize()
            })
        },
        destroy: function() {
            this.autoPlace && T.removeChild(this.domElement)
        },
        addFolder: function(e) {
            if (void 0 !== this.__folders[e])
                throw new Error('You already have a folder in this GUI by the name "' + e + '"');
            var t = {
                name: e,
                parent: this
            };
            t.autoPlace = this.autoPlace,
            this.load && this.load.folders && this.load.folders[e] && (t.closed = this.load.folders[e].closed,
            t.load = this.load.folders[e]);
            var n = new H(t);
            this.__folders[e] = n;
            var o = m(this, n.domElement);
            return p.addClass(o, "folder"),
            n
        },
        open: function() {
            this.closed = !1
        },
        close: function() {
            this.closed = !0
        },
        onResize: function() {
            var e = this.getRoot();
            if (e.scrollable) {
                var t = p.getOffset(e.__ul).top
                  , n = 0;
                _.each(e.__ul.childNodes, function(t) {
                    e.autoPlace && t === e.__save_row || (n += p.getHeight(t))
                }),
                window.innerHeight - t - R < n ? (p.addClass(e.domElement, H.CLASS_TOO_TALL),
                e.__ul.style.height = window.innerHeight - t - R + "px") : (p.removeClass(e.domElement, H.CLASS_TOO_TALL),
                e.__ul.style.height = "auto")
            }
            e.__resize_handle && _.defer(function() {
                e.__resize_handle.style.height = e.__ul.offsetHeight + "px"
            }),
            e.__closeButton && (e.__closeButton.style.width = e.width + "px")
        },
        remember: function() {
            if (_.isUndefined(O) && (O = new h,
            O.domElement.innerHTML = t),
            this.parent)
                throw new Error("You can only call remember on a top level GUI.");
            var e = this;
            _.each(Array.prototype.slice.call(arguments), function(t) {
                0 == e.__rememberedObjects.length && y(e),
                -1 == e.__rememberedObjects.indexOf(t) && e.__rememberedObjects.push(t)
            }),
            this.autoPlace && w(this, this.width)
        },
        getRoot: function() {
            for (var e = this; e.parent; )
                e = e.parent;
            return e
        },
        getSaveObject: function() {
            var e = this.load;
            return e.closed = this.closed,
            this.__rememberedObjects.length > 0 && (e.preset = this.preset,
            e.remembered || (e.remembered = {}),
            e.remembered[this.preset] = C(this)),
            e.folders = {},
            _.each(this.__folders, function(t, n) {
                e.folders[n] = t.getSaveObject()
            }),
            e
        },
        save: function() {
            this.load.remembered || (this.load.remembered = {}),
            this.load.remembered[this.preset] = C(this),
            k(this, !1),
            this.saveToLocalStorageIfPossible()
        },
        saveAs: function(e) {
            this.load.remembered || (this.load.remembered = {},
            this.load.remembered[B] = C(this, !0)),
            this.load.remembered[e] = C(this),
            this.preset = e,
            E(this, e, !0),
            this.saveToLocalStorageIfPossible()
        },
        revert: function(e) {
            _.each(this.__controllers, function(t) {
                this.getRoot().load.remembered ? b(e || this.getRoot(), t) : t.setValue(t.initialValue)
            }, this),
            _.each(this.__folders, function(e) {
                e.revert(e)
            }),
            e || k(this.getRoot(), !1)
        },
        listen: function(e) {
            var t = 0 == this.__listening.length;
            this.__listening.push(e),
            t && S(this.__listening)
        }
    }),
    H
}(dat.utils.css, '<div id="dg-save" class="dg dialogue">\n\n  Here\'s the new load parameter for your <code>GUI</code>\'s constructor:\n\n  <textarea id="dg-new-constructor"></textarea>\n\n  <div id="dg-save-locally">\n\n    <input id="dg-local-storage" type="checkbox"/> Automatically save\n    values to <code>localStorage</code> on exit.\n\n    <div id="dg-local-explain">The values saved to <code>localStorage</code> will\n      override those passed to <code>dat.GUI</code>\'s constructor. This makes it\n      easier to work incrementally, but <code>localStorage</code> is fragile,\n      and your friends may not see the same values you do.\n      \n    </div>\n    \n  </div>\n\n</div>', ".dg {\n  /** Clear list styles */\n  /* Auto-place container */\n  /* Auto-placed GUI's */\n  /* Line items that don't contain folders. */\n  /** Folder names */\n  /** Hides closed items */\n  /** Controller row */\n  /** Name-half (left) */\n  /** Controller-half (right) */\n  /** Controller placement */\n  /** Shorter number boxes when slider is present. */\n  /** Ensure the entire boolean and function row shows a hand */ }\n  .dg ul {\n    list-style: none;\n    margin: 0;\n    padding: 0;\n    width: 100%;\n    clear: both; }\n  .dg.ac {\n    position: fixed;\n    top: 0;\n    left: 0;\n    right: 0;\n    height: 0;\n    z-index: 0; }\n  .dg:not(.ac) .main {\n    /** Exclude mains in ac so that we don't hide close button */\n    overflow: hidden; }\n  .dg.main {\n    -webkit-transition: opacity 0.1s linear;\n    -o-transition: opacity 0.1s linear;\n    -moz-transition: opacity 0.1s linear;\n    transition: opacity 0.1s linear; }\n    .dg.main.taller-than-window {\n      overflow-y: auto; }\n      .dg.main.taller-than-window .close-button {\n        opacity: 1;\n        /* TODO, these are style notes */\n        margin-top: -1px;\n        border-top: 1px solid #2c2c2c; }\n    .dg.main ul.closed .close-button {\n      opacity: 1 !important; }\n    .dg.main:hover .close-button,\n    .dg.main .close-button.drag {\n      opacity: 1; }\n    .dg.main .close-button {\n      /*opacity: 0;*/\n      -webkit-transition: opacity 0.1s linear;\n      -o-transition: opacity 0.1s linear;\n      -moz-transition: opacity 0.1s linear;\n      transition: opacity 0.1s linear;\n      border: 0;\n      position: absolute;\n      line-height: 19px;\n      height: 20px;\n      /* TODO, these are style notes */\n      cursor: pointer;\n      text-align: center;\n      background-color: #000; }\n      .dg.main .close-button:hover {\n        background-color: #111; }\n  .dg.a {\n    float: right;\n    margin-right: 15px;\n    overflow-x: hidden; }\n    .dg.a.has-save > ul {\n      margin-top: 27px; }\n      .dg.a.has-save > ul.closed {\n        margin-top: 0; }\n    .dg.a .save-row {\n      position: fixed;\n      top: 0;\n      z-index: 1002; }\n  .dg li {\n    -webkit-transition: height 0.1s ease-out;\n    -o-transition: height 0.1s ease-out;\n    -moz-transition: height 0.1s ease-out;\n    transition: height 0.1s ease-out; }\n  .dg li:not(.folder) {\n    cursor: auto;\n    height: 27px;\n    line-height: 27px;\n    overflow: hidden;\n    padding: 0 4px 0 5px; }\n  .dg li.folder {\n    padding: 0;\n    border-left: 4px solid rgba(0, 0, 0, 0); }\n  .dg li.title {\n    cursor: pointer;\n    margin-left: -4px; }\n  .dg .closed li:not(.title),\n  .dg .closed ul li,\n  .dg .closed ul li > * {\n    height: 0;\n    overflow: hidden;\n    border: 0; }\n  .dg .cr {\n    clear: both;\n    padding-left: 3px;\n    height: 27px; }\n  .dg .property-name {\n    cursor: default;\n    float: left;\n    clear: left;\n    width: 40%;\n    overflow: hidden;\n    text-overflow: ellipsis; }\n  .dg .c {\n    float: left;\n    width: 60%; }\n  .dg .c input[type=text] {\n    border: 0;\n    margin-top: 4px;\n    padding: 3px;\n    width: 100%;\n    float: right; }\n  .dg .has-slider input[type=text] {\n    width: 30%;\n    /*display: none;*/\n    margin-left: 0; }\n  .dg .slider {\n    float: left;\n    width: 66%;\n    margin-left: -5px;\n    margin-right: 0;\n    height: 19px;\n    margin-top: 4px; }\n  .dg .slider-fg {\n    height: 100%; }\n  .dg .c input[type=checkbox] {\n    margin-top: 9px; }\n  .dg .c select {\n    margin-top: 5px; }\n  .dg .cr.function,\n  .dg .cr.function .property-name,\n  .dg .cr.function *,\n  .dg .cr.boolean,\n  .dg .cr.boolean * {\n    cursor: pointer; }\n  .dg .selector {\n    display: none;\n    position: absolute;\n    margin-left: -9px;\n    margin-top: 23px;\n    z-index: 10; }\n  .dg .c:hover .selector,\n  .dg .selector.drag {\n    display: block; }\n  .dg li.save-row {\n    padding: 0; }\n    .dg li.save-row .button {\n      display: inline-block;\n      padding: 0px 6px; }\n  .dg.dialogue {\n    background-color: #222;\n    width: 460px;\n    padding: 15px;\n    font-size: 13px;\n    line-height: 15px; }\n\n/* TODO Separate style and structure */\n#dg-new-constructor {\n  padding: 10px;\n  color: #222;\n  font-family: Monaco, monospace;\n  font-size: 10px;\n  border: 0;\n  resize: none;\n  box-shadow: inset 1px 1px 1px #888;\n  word-wrap: break-word;\n  margin: 12px 0;\n  display: block;\n  width: 440px;\n  overflow-y: scroll;\n  height: 100px;\n  position: relative; }\n\n#dg-local-explain {\n  display: none;\n  font-size: 11px;\n  line-height: 17px;\n  border-radius: 3px;\n  background-color: #333;\n  padding: 8px;\n  margin-top: 10px; }\n  #dg-local-explain code {\n    font-size: 10px; }\n\n#dat-gui-save-locally {\n  display: none; }\n\n/** Main type */\n.dg {\n  color: #eee;\n  font: 11px 'Lucida Grande', sans-serif;\n  text-shadow: 0 -1px 0 #111;\n  /** Auto place */\n  /* Controller row, <li> */\n  /** Controllers */ }\n  .dg.main {\n    /** Scrollbar */ }\n    .dg.main::-webkit-scrollbar {\n      width: 5px;\n      background: #1a1a1a; }\n    .dg.main::-webkit-scrollbar-corner {\n      height: 0;\n      display: none; }\n    .dg.main::-webkit-scrollbar-thumb {\n      border-radius: 5px;\n      background: #676767; }\n  .dg li:not(.folder) {\n    background: #1a1a1a;\n    border-bottom: 1px solid #2c2c2c; }\n  .dg li.save-row {\n    line-height: 25px;\n    background: #dad5cb;\n    border: 0; }\n    .dg li.save-row select {\n      margin-left: 5px;\n      width: 108px; }\n    .dg li.save-row .button {\n      margin-left: 5px;\n      margin-top: 1px;\n      border-radius: 2px;\n      font-size: 9px;\n      line-height: 7px;\n      padding: 4px 4px 5px 4px;\n      background: #c5bdad;\n      color: #fff;\n      text-shadow: 0 1px 0 #b0a58f;\n      box-shadow: 0 -1px 0 #b0a58f;\n      cursor: pointer; }\n      .dg li.save-row .button.gears {\n        background: #c5bdad url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAsAAAANCAYAAAB/9ZQ7AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAQJJREFUeNpiYKAU/P//PwGIC/ApCABiBSAW+I8AClAcgKxQ4T9hoMAEUrxx2QSGN6+egDX+/vWT4e7N82AMYoPAx/evwWoYoSYbACX2s7KxCxzcsezDh3evFoDEBYTEEqycggWAzA9AuUSQQgeYPa9fPv6/YWm/Acx5IPb7ty/fw+QZblw67vDs8R0YHyQhgObx+yAJkBqmG5dPPDh1aPOGR/eugW0G4vlIoTIfyFcA+QekhhHJhPdQxbiAIguMBTQZrPD7108M6roWYDFQiIAAv6Aow/1bFwXgis+f2LUAynwoIaNcz8XNx3Dl7MEJUDGQpx9gtQ8YCueB+D26OECAAQDadt7e46D42QAAAABJRU5ErkJggg==) 2px 1px no-repeat;\n        height: 7px;\n        width: 8px; }\n      .dg li.save-row .button:hover {\n        background-color: #bab19e;\n        box-shadow: 0 -1px 0 #b0a58f; }\n  .dg li.folder {\n    border-bottom: 0; }\n  .dg li.title {\n    padding-left: 16px;\n    background: black url(data:image/gif;base64,R0lGODlhBQAFAJEAAP////Pz8////////yH5BAEAAAIALAAAAAAFAAUAAAIIlI+hKgFxoCgAOw==) 6px 10px no-repeat;\n    cursor: pointer;\n    border-bottom: 1px solid rgba(255, 255, 255, 0.2); }\n  .dg .closed li.title {\n    background-image: url(data:image/gif;base64,R0lGODlhBQAFAJEAAP////Pz8////////yH5BAEAAAIALAAAAAAFAAUAAAIIlGIWqMCbWAEAOw==); }\n  .dg .cr.boolean {\n    border-left: 3px solid #806787; }\n  .dg .cr.function {\n    border-left: 3px solid #e61d5f; }\n  .dg .cr.number {\n    border-left: 3px solid #2fa1d6; }\n    .dg .cr.number input[type=text] {\n      color: #2fa1d6; }\n  .dg .cr.string {\n    border-left: 3px solid #1ed36f; }\n    .dg .cr.string input[type=text] {\n      color: #1ed36f; }\n  .dg .cr.function:hover, .dg .cr.boolean:hover {\n    background: #111; }\n  .dg .c input[type=text] {\n    background: #303030;\n    outline: none; }\n    .dg .c input[type=text]:hover {\n      background: #3c3c3c; }\n    .dg .c input[type=text]:focus {\n      background: #494949;\n      color: #fff; }\n  .dg .c .slider {\n    background: #303030;\n    cursor: ew-resize; }\n  .dg .c .slider-fg {\n    background: #2fa1d6; }\n  .dg .c .slider:hover {\n    background: #3c3c3c; }\n    .dg .c .slider:hover .slider-fg {\n      background: #44abda; }\n", dat.controllers.factory = function(e, t, n, o, i, r, s) {
    return function(a, l) {
        var d = a[l];
        return s.isArray(arguments[2]) || s.isObject(arguments[2]) ? new e(a,l,arguments[2]) : s.isNumber(d) ? s.isNumber(arguments[2]) && s.isNumber(arguments[3]) ? new n(a,l,arguments[2],arguments[3]) : new t(a,l,{
            min: arguments[2],
            max: arguments[3]
        }) : s.isString(d) ? new o(a,l) : s.isFunction(d) ? new i(a,l,"") : s.isBoolean(d) ? new r(a,l) : void 0
    }
}(dat.controllers.OptionController, dat.controllers.NumberControllerBox, dat.controllers.NumberControllerSlider, dat.controllers.StringController = function(e, t, n) {
    var o = function(e, n) {
        function i() {
            s.setValue(s.__input.value)
        }
        function r() {
            s.__onFinishChange && s.__onFinishChange.call(s, s.getValue())
        }
        o.superclass.call(this, e, n);
        var s = this;
        this.__input = document.createElement("input"),
        this.__input.setAttribute("type", "text"),
        t.bind(this.__input, "keyup", i),
        t.bind(this.__input, "change", i),
        t.bind(this.__input, "blur", r),
        t.bind(this.__input, "keydown", function(e) {
            13 === e.keyCode && this.blur()
        }),
        this.updateDisplay(),
        this.domElement.appendChild(this.__input)
    };
    return o.superclass = e,
    n.extend(o.prototype, e.prototype, {
        updateDisplay: function() {
            return t.isActive(this.__input) || (this.__input.value = this.getValue()),
            o.superclass.prototype.updateDisplay.call(this)
        }
    }),
    o
}(dat.controllers.Controller, dat.dom.dom, dat.utils.common), dat.controllers.FunctionController, dat.controllers.BooleanController, dat.utils.common), dat.controllers.Controller, dat.controllers.BooleanController, dat.controllers.FunctionController, dat.controllers.NumberControllerBox, dat.controllers.NumberControllerSlider, dat.controllers.OptionController, dat.controllers.ColorController = function(e, t, n, o, i) {
    function r(e, t, n, o) {
        e.style.background = "",
        i.each(l, function(i) {
            e.style.cssText += "background: " + i + "linear-gradient(" + t + ", " + n + " 0%, " + o + " 100%); "
        })
    }
    function s(e) {
        e.style.background = "",
        e.style.cssText += "background: -moz-linear-gradient(top,  #ff0000 0%, #ff00ff 17%, #0000ff 34%, #00ffff 50%, #00ff00 67%, #ffff00 84%, #ff0000 100%);",
        e.style.cssText += "background: -webkit-linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);",
        e.style.cssText += "background: -o-linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);",
        e.style.cssText += "background: -ms-linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);",
        e.style.cssText += "background: linear-gradient(top,  #ff0000 0%,#ff00ff 17%,#0000ff 34%,#00ffff 50%,#00ff00 67%,#ffff00 84%,#ff0000 100%);"
    }
    var a = function(e, l) {
        function d(e) {
            p(e),
            t.bind(window, "mousemove", p),
            t.bind(window, "mouseup", c)
        }
        function c() {
            t.unbind(window, "mousemove", p),
            t.unbind(window, "mouseup", c)
        }
        function u() {
            var e = o(this.value);
            e !== !1 ? (f.__color.__state = e,
            f.setValue(f.__color.toOriginal())) : this.value = f.__color.toString()
        }
        function h() {
            t.unbind(window, "mousemove", _),
            t.unbind(window, "mouseup", h)
        }
        function p(e) {
            e.preventDefault();
            var n = t.getWidth(f.__saturation_field)
              , o = t.getOffset(f.__saturation_field)
              , i = (e.clientX - o.left + document.body.scrollLeft) / n
              , r = 1 - (e.clientY - o.top + document.body.scrollTop) / n;
            return r > 1 ? r = 1 : 0 > r && (r = 0),
            i > 1 ? i = 1 : 0 > i && (i = 0),
            f.__color.v = r,
            f.__color.s = i,
            f.setValue(f.__color.toOriginal()),
            !1
        }
        function _(e) {
            e.preventDefault();
            var n = t.getHeight(f.__hue_field)
              , o = t.getOffset(f.__hue_field)
              , i = 1 - (e.clientY - o.top + document.body.scrollTop) / n;
            return i > 1 ? i = 1 : 0 > i && (i = 0),
            f.__color.h = 360 * i,
            f.setValue(f.__color.toOriginal()),
            !1
        }
        a.superclass.call(this, e, l),
        this.__color = new n(this.getValue()),
        this.__temp = new n(0);
        var f = this;
        this.domElement = document.createElement("div"),
        t.makeSelectable(this.domElement, !1),
        this.__selector = document.createElement("div"),
        this.__selector.className = "selector",
        this.__saturation_field = document.createElement("div"),
        this.__saturation_field.className = "saturation-field",
        this.__field_knob = document.createElement("div"),
        this.__field_knob.className = "field-knob",
        this.__field_knob_border = "2px solid ",
        this.__hue_knob = document.createElement("div"),
        this.__hue_knob.className = "hue-knob",
        this.__hue_field = document.createElement("div"),
        this.__hue_field.className = "hue-field",
        this.__input = document.createElement("input"),
        this.__input.type = "text",
        this.__input_textShadow = "0 1px 1px ",
        t.bind(this.__input, "keydown", function(e) {
            13 === e.keyCode && u.call(this)
        }),
        t.bind(this.__input, "blur", u),
        t.bind(this.__selector, "mousedown", function() {
            t.addClass(this, "drag").bind(window, "mouseup", function() {
                t.removeClass(f.__selector, "drag")
            })
        });
        var m = document.createElement("div");
        i.extend(this.__selector.style, {
            width: "122px",
            height: "102px",
            padding: "3px",
            backgroundColor: "#222",
            boxShadow: "0px 1px 3px rgba(0,0,0,0.3)"
        }),
        i.extend(this.__field_knob.style, {
            position: "absolute",
            width: "12px",
            height: "12px",
            border: this.__field_knob_border + (this.__color.v < .5 ? "#fff" : "#000"),
            boxShadow: "0px 1px 3px rgba(0,0,0,0.5)",
            borderRadius: "12px",
            zIndex: 1
        }),
        i.extend(this.__hue_knob.style, {
            position: "absolute",
            width: "15px",
            height: "2px",
            borderRight: "4px solid #fff",
            zIndex: 1
        }),
        i.extend(this.__saturation_field.style, {
            width: "100px",
            height: "100px",
            border: "1px solid #555",
            marginRight: "3px",
            display: "inline-block",
            cursor: "pointer"
        }),
        i.extend(m.style, {
            width: "100%",
            height: "100%",
            background: "none"
        }),
        r(m, "top", "rgba(0,0,0,0)", "#000"),
        i.extend(this.__hue_field.style, {
            width: "15px",
            height: "100px",
            display: "inline-block",
            border: "1px solid #555",
            cursor: "ns-resize"
        }),
        s(this.__hue_field),
        i.extend(this.__input.style, {
            outline: "none",
            textAlign: "center",
            color: "#fff",
            border: 0,
            fontWeight: "bold",
            textShadow: this.__input_textShadow + "rgba(0,0,0,0.7)"
        }),
        t.bind(this.__saturation_field, "mousedown", d),
        t.bind(this.__field_knob, "mousedown", d),
        t.bind(this.__hue_field, "mousedown", function(e) {
            _(e),
            t.bind(window, "mousemove", _),
            t.bind(window, "mouseup", h)
        }),
        this.__saturation_field.appendChild(m),
        this.__selector.appendChild(this.__field_knob),
        this.__selector.appendChild(this.__saturation_field),
        this.__selector.appendChild(this.__hue_field),
        this.__hue_field.appendChild(this.__hue_knob),
        this.domElement.appendChild(this.__input),
        this.domElement.appendChild(this.__selector),
        this.updateDisplay()
    };
    a.superclass = e,
    i.extend(a.prototype, e.prototype, {
        updateDisplay: function() {
            var e = o(this.getValue());
            if (e !== !1) {
                var t = !1;
                i.each(n.COMPONENTS, function(n) {
                    return i.isUndefined(e[n]) || i.isUndefined(this.__color.__state[n]) || e[n] === this.__color.__state[n] ? void 0 : (t = !0,
                    {})
                }, this),
                t && i.extend(this.__color.__state, e)
            }
            i.extend(this.__temp.__state, this.__color.__state),
            this.__temp.a = 1;
            var s = this.__color.v < .5 || this.__color.s > .5 ? 255 : 0
              , a = 255 - s;
            i.extend(this.__field_knob.style, {
                marginLeft: 100 * this.__color.s - 7 + "px",
                marginTop: 100 * (1 - this.__color.v) - 7 + "px",
                backgroundColor: this.__temp.toString(),
                border: this.__field_knob_border + "rgb(" + s + "," + s + "," + s + ")"
            }),
            this.__hue_knob.style.marginTop = 100 * (1 - this.__color.h / 360) + "px",
            this.__temp.s = 1,
            this.__temp.v = 1,
            r(this.__saturation_field, "left", "#fff", this.__temp.toString()),
            i.extend(this.__input.style, {
                backgroundColor: this.__input.value = this.__color.toString(),
                color: "rgb(" + s + "," + s + "," + s + ")",
                textShadow: this.__input_textShadow + "rgba(" + a + "," + a + "," + a + ",.7)"
            })
        }
    });
    var l = ["-moz-", "-o-", "-webkit-", "-ms-", ""];
    return a
}(dat.controllers.Controller, dat.dom.dom, dat.color.Color = function(e, t, n, o) {
    function i(e, t, n) {
        Object.defineProperty(e, t, {
            get: function() {
                return "RGB" === this.__state.space ? this.__state[t] : (s(this, t, n),
                this.__state[t])
            },
            set: function(e) {
                "RGB" !== this.__state.space && (s(this, t, n),
                this.__state.space = "RGB"),
                this.__state[t] = e
            }
        })
    }
    function r(e, t) {
        Object.defineProperty(e, t, {
            get: function() {
                return "HSV" === this.__state.space ? this.__state[t] : (a(this),
                this.__state[t])
            },
            set: function(e) {
                "HSV" !== this.__state.space && (a(this),
                this.__state.space = "HSV"),
                this.__state[t] = e
            }
        })
    }
    function s(e, n, i) {
        if ("HEX" === e.__state.space)
            e.__state[n] = t.component_from_hex(e.__state.hex, i);
        else {
            if ("HSV" !== e.__state.space)
                throw "Corrupted color state";
            o.extend(e.__state, t.hsv_to_rgb(e.__state.h, e.__state.s, e.__state.v))
        }
    }
    function a(e) {
        var n = t.rgb_to_hsv(e.r, e.g, e.b);
        o.extend(e.__state, {
            s: n.s,
            v: n.v
        }),
        o.isNaN(n.h) ? o.isUndefined(e.__state.h) && (e.__state.h = 0) : e.__state.h = n.h
    }
    var l = function() {
        if (this.__state = e.apply(this, arguments),
        this.__state === !1)
            throw "Failed to interpret color arguments";
        this.__state.a = this.__state.a || 1
    };
    return l.COMPONENTS = ["r", "g", "b", "h", "s", "v", "hex", "a"],
    o.extend(l.prototype, {
        toString: function() {
            return n(this)
        },
        toOriginal: function() {
            return this.__state.conversion.write(this)
        }
    }),
    i(l.prototype, "r", 2),
    i(l.prototype, "g", 1),
    i(l.prototype, "b", 0),
    r(l.prototype, "h"),
    r(l.prototype, "s"),
    r(l.prototype, "v"),
    Object.defineProperty(l.prototype, "a", {
        get: function() {
            return this.__state.a
        },
        set: function(e) {
            this.__state.a = e
        }
    }),
    Object.defineProperty(l.prototype, "hex", {
        get: function() {
            return "HEX" !== !this.__state.space && (this.__state.hex = t.rgb_to_hex(this.r, this.g, this.b)),
            this.__state.hex
        },
        set: function(e) {
            this.__state.space = "HEX",
            this.__state.hex = e
        }
    }),
    l
}(dat.color.interpret, dat.color.math = function() {
    var e;
    return {
        hsv_to_rgb: function(e, t, n) {
            var o = Math.floor(e / 60) % 6
              , i = e / 60 - Math.floor(e / 60)
              , r = n * (1 - t)
              , s = n * (1 - i * t)
              , a = n * (1 - (1 - i) * t)
              , l = [[n, a, r], [s, n, r], [r, n, a], [r, s, n], [a, r, n], [n, r, s]][o];
            return {
                r: 255 * l[0],
                g: 255 * l[1],
                b: 255 * l[2]
            }
        },
        rgb_to_hsv: function(e, t, n) {
            var o, i, r = Math.min(e, t, n), s = Math.max(e, t, n), a = s - r;
            return 0 == s ? {
                h: 0 / 0,
                s: 0,
                v: 0
            } : (i = a / s,
            o = e == s ? (t - n) / a : t == s ? 2 + (n - e) / a : 4 + (e - t) / a,
            o /= 6,
            0 > o && (o += 1),
            {
                h: 360 * o,
                s: i,
                v: s / 255
            })
        },
        rgb_to_hex: function(e, t, n) {
            var o = this.hex_with_component(0, 2, e);
            return o = this.hex_with_component(o, 1, t),
            o = this.hex_with_component(o, 0, n)
        },
        component_from_hex: function(e, t) {
            return e >> 8 * t & 255
        },
        hex_with_component: function(t, n, o) {
            return o << (e = 8 * n) | t & ~(255 << e)
        }
    }
}(), dat.color.toString, dat.utils.common), dat.color.interpret, dat.utils.common), dat.utils.requestAnimationFrame = function() {
    return window.webkitRequestAnimationFrame || window.mozRequestAnimationFrame || window.oRequestAnimationFrame || window.msRequestAnimationFrame || function(e) {
        window.setTimeout(e, 1e3 / 60)
    }
}(), dat.dom.CenteredDiv = function(e, t) {
    var n = function() {
        this.backgroundElement = document.createElement("div"),
        t.extend(this.backgroundElement.style, {
            backgroundColor: "rgba(0,0,0,0.8)",
            top: 0,
            left: 0,
            display: "none",
            zIndex: "1000",
            opacity: 0,
            WebkitTransition: "opacity 0.2s linear",
            transition: "opacity 0.2s linear"
        }),
        e.makeFullscreen(this.backgroundElement),
        this.backgroundElement.style.position = "fixed",
        this.domElement = document.createElement("div"),
        t.extend(this.domElement.style, {
            position: "fixed",
            display: "none",
            zIndex: "1001",
            opacity: 0,
            WebkitTransition: "-webkit-transform 0.2s ease-out, opacity 0.2s linear",
            transition: "transform 0.2s ease-out, opacity 0.2s linear"
        }),
        document.body.appendChild(this.backgroundElement),
        document.body.appendChild(this.domElement);
        var n = this;
        e.bind(this.backgroundElement, "click", function() {
            n.hide()
        })
    };
    return n.prototype.show = function() {
        var e = this;
        this.backgroundElement.style.display = "block",
        this.domElement.style.display = "block",
        this.domElement.style.opacity = 0,
        this.domElement.style.webkitTransform = "scale(1.1)",
        this.layout(),
        t.defer(function() {
            e.backgroundElement.style.opacity = 1,
            e.domElement.style.opacity = 1,
            e.domElement.style.webkitTransform = "scale(1)"
        })
    }
    ,
    n.prototype.hide = function() {
        var t = this
          , n = function() {
            t.domElement.style.display = "none",
            t.backgroundElement.style.display = "none",
            e.unbind(t.domElement, "webkitTransitionEnd", n),
            e.unbind(t.domElement, "transitionend", n),
            e.unbind(t.domElement, "oTransitionEnd", n)
        };
        e.bind(this.domElement, "webkitTransitionEnd", n),
        e.bind(this.domElement, "transitionend", n),
        e.bind(this.domElement, "oTransitionEnd", n),
        this.backgroundElement.style.opacity = 0,
        this.domElement.style.opacity = 0,
        this.domElement.style.webkitTransform = "scale(1.1)"
    }
    ,
    n.prototype.layout = function() {
        this.domElement.style.left = window.innerWidth / 2 - e.getWidth(this.domElement) / 2 + "px",
        this.domElement.style.top = window.innerHeight / 2 - e.getHeight(this.domElement) / 2 + "px"
    }
    ,
    n
}(dat.dom.dom, dat.utils.common), dat.dom.dom, dat.utils.common);
