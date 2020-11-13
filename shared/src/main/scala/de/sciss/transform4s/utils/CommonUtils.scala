/* ***** BEGIN LICENSE BLOCK *****
 * JTransforms
 * Copyright (c) 2007 onward, Piotr Wendykier
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ***** END LICENSE BLOCK ***** */

package de.sciss.transform4s.utils

import scala.concurrent.Future
//import org.apache.commons.math3.util.FastMath._
import ConcurrencyUtils.executionContext

import Math.{floor, cos, sin, log, pow}

/**
 * Static utility methods.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
object CommonUtils {
  private var THREADS_BEGIN_N_1D_FFT_2THREADS: Int = 8192
  private var THREADS_BEGIN_N_1D_FFT_4THREADS: Int = 65536
  private var THREADS_BEGIN_N_2D             : Int = 65536
  private var THREADS_BEGIN_N_3D             : Int = 65536

  /**
   * Causes the currently executing thread to sleep (temporarily cease
   * execution) for the specified number of milliseconds.
   *
   * @param millis the length of time to sleep in milliseconds
   */
  def sleep(millis: Long): Unit = {
    try Thread.sleep(millis)
    catch {
      case e: InterruptedException =>
        e.printStackTrace()
    }
  }

  /**
   * Returns the minimal size of 1D data for which two threads are used.
   *
   * @return the minimal size of 1D data for which two threads are used
   */
  def threadsBeginN_1D_FFT_2Threads: Int = THREADS_BEGIN_N_1D_FFT_2THREADS

  /**
   * Returns the minimal size of 1D data for which four threads are used.
   *
   * @return the minimal size of 1D data for which four threads are used
   */
  def threadsBeginN_1D_FFT_4Threads: Int = THREADS_BEGIN_N_1D_FFT_4THREADS

  /**
   * Returns the minimal size of 2D data for which threads are used.
   *
   * @return the minimal size of 2D data for which threads are used
   */
  def threadsBeginN_2D: Int = THREADS_BEGIN_N_2D

  /**
   * Returns the minimal size of 3D data for which threads are used.
   *
   * @return the minimal size of 3D data for which threads are used
   */
  def threadsBeginN_3D: Int = THREADS_BEGIN_N_3D

  /**
   * Sets the minimal size of 1D data for which two threads are used.
   *
   * @param n the minimal size of 1D data for which two threads are used
   */
  def threadsBeginN_1D_FFT_2Threads_=(n: Int): Unit =
    if (n < 1024) THREADS_BEGIN_N_1D_FFT_2THREADS = 1024
    else          THREADS_BEGIN_N_1D_FFT_2THREADS = n

  /**
   * Sets the minimal size of 1D data for which four threads are used.
   *
   * @param n the minimal size of 1D data for which four threads are used
   */
  def threadsBeginN_1D_FFT_4Threads_=(n: Int): Unit =
    if (n < 1024) THREADS_BEGIN_N_1D_FFT_4THREADS = 1024
    else          THREADS_BEGIN_N_1D_FFT_4THREADS = n

  /**
   * Sets the minimal size of 2D data for which threads are used.
   *
   * @param n the minimal size of 2D data for which threads are used
   */
  def threadsBeginN_2D_=(n: Int): Unit =
    if (n < 4096) THREADS_BEGIN_N_2D = 4096
    else          THREADS_BEGIN_N_2D = n

  /**
   * Sets the minimal size of 3D data for which threads are used.
   *
   * @param n the minimal size of 3D data for which threads are used
   */
  def threadsBeginN_3D_=(n: Int): Unit =
    THREADS_BEGIN_N_3D = n

  /**
   * Resets the minimal size of 1D data for which two and four threads are
   * used.
   */
  def resetThreadsBeginN_FFT(): Unit = {
    THREADS_BEGIN_N_1D_FFT_2THREADS = 8192
    THREADS_BEGIN_N_1D_FFT_4THREADS = 65536
  }

  /**
   * Resets the minimal size of 2D and 3D data for which threads are used.
   */
  def resetThreadsBeginN(): Unit = {
    THREADS_BEGIN_N_2D = 65536
    THREADS_BEGIN_N_3D = 65536
  }

  /**
   * Returns the closest power-of-two number greater than or equal to x.
   *
   * @param x input value
   * @return the closest power-of-two number greater than or equal to x
   */
  def nextPow2(x: Int): Int = {
    if (x < 1) throw new IllegalArgumentException("x must be greater or equal 1")
    if ((x & (x - 1)) == 0) return x // x is already a power-of-two number
    var xi = x
    xi |= (xi >>>  1)
    xi |= (xi >>>  2)
    xi |= (xi >>>  4)
    xi |= (xi >>>  8)
    xi |= (xi >>> 16)
    xi + 1
  }

  def nextPow2(x: Long): Long = {
    if (x < 1) throw new IllegalArgumentException("x must be greater or equal 1")
    if ((x & (x - 1L)) == 0) return x
    var xi = x
    xi |= (xi >>>  1L)
    xi |= (xi >>>  2L)
    xi |= (xi >>>  4L)
    xi |= (xi >>>  8L)
    xi |= (xi >>> 16L)
    xi |= (xi >>> 32L)
    xi + 1L
  }

  /**
   * Returns the closest power-of-two number less than or equal to x.
   *
   * @param x input value
   * @return the closest power-of-two number less then or equal to x
   */
  def prevPow2(x: Int): Int = {
    if (x < 1) throw new IllegalArgumentException("x must be greater or equal 1")
    pow(2, floor(log(x) / log(2))).toInt
  }

  def prevPow2(x: Long): Long = {
    if (x < 1) throw new IllegalArgumentException("x must be greater or equal 1")
    pow(2, floor(log(x.toDouble) / log(2))).toLong
  }

  /**
   * Checks if x is a power-of-two number.
   *
   * @param x input value
   * @return true if x is a power-of-two number
   */
  def isPowerOf2(x: Int): Boolean = if (x <= 0) false
  else (x & (x - 1)) == 0

  def isPowerOf2(x: Long): Boolean = if (x <= 0) false
  else (x & (x - 1L)) == 0

  def getReminder(n: Long, factors: Array[Int]): Long = {
    var reminder = n
    if (n <= 0) throw new IllegalArgumentException("n must be positive integer")
    var i = 0
    while ( {
      i < factors.length && reminder != 1L
    }) {
      val factor = factors(i)
      while ( {
        (reminder % factor) == 0
      }) reminder /= factor

      i += 1
    }
    reminder
  }

  def makeipt(nw: Int, ip: Array[Int]): Unit = {
    var j = 0
    var l = 0
    var m = 0
    var m2 = 0
    var p = 0
    var q = 0
    ip(2) = 0
    ip(3) = 16
    m = 2
    l = nw
    while ( {
      l > 32
    }) {
      m2 = m << 1
      q = m2 << 3
      j = m
      while ( {
        j < m2
      }) {
        p = ip(j) << 2
        ip(m + j) = p
        ip(m2 + j) = p + q

        j += 1
      }
      m = m2

      l >>= 2
    }
  }

  def makewt(nw: Int, ip: Array[Int], w: Array[Double]): Unit = {
    var j = 0
    var nwh = 0
    var nw0 = 0
    var nw1 = 0
    var delta = .0
    var wn4r = .0
    var wk1r = .0
    var wk1i = .0
    var wk3r = .0
    var wk3i = .0
    var delta2 = .0
    var deltaj = .0
    var deltaj3 = .0
    ip(0) = nw
    ip(1) = 1
    if (nw > 2) {
      nwh = nw >> 1
      delta = 0.785398163397448278999490867136046290 / nwh
      delta2 = delta * 2
      wn4r = cos(delta * nwh)
      w(0) = 1
      w(1) = wn4r
      if (nwh == 4) {
        w(2) = cos(delta2)
        w(3) = sin(delta2)
      }
      else if (nwh > 4) {
        CommonUtils.makeipt(nw, ip)
        w(2) = 0.5 / cos(delta2)
        w(3) = 0.5 / cos(delta * 6)
        j = 4
        while ( {
          j < nwh
        }) {
          deltaj = delta * j
          deltaj3 = 3 * deltaj
          w(j) = cos(deltaj)
          w(j + 1) = sin(deltaj)
          w(j + 2) = cos(deltaj3)
          w(j + 3) = -sin(deltaj3)

          j += 4
        }
      }
      nw0 = 0
      while ( {
        nwh > 2
      }) {
        nw1 = nw0 + nwh
        nwh >>= 1
        w(nw1) = 1
        w(nw1 + 1) = wn4r
        if (nwh == 4) {
          wk1r = w(nw0 + 4)
          wk1i = w(nw0 + 5)
          w(nw1 + 2) = wk1r
          w(nw1 + 3) = wk1i
        }
        else if (nwh > 4) {
          wk1r = w(nw0 + 4)
          wk3r = w(nw0 + 6)
          w(nw1 + 2) = 0.5 / wk1r
          w(nw1 + 3) = 0.5 / wk3r
          j = 4
          while ( {
            j < nwh
          }) {
            val idx1 = nw0 + 2 * j
            val idx2 = nw1 + j
            wk1r = w(idx1)
            wk1i = w(idx1 + 1)
            wk3r = w(idx1 + 2)
            wk3i = w(idx1 + 3)
            w(idx2) = wk1r
            w(idx2 + 1) = wk1i
            w(idx2 + 2) = wk3r
            w(idx2 + 3) = wk3i

            j += 4
          }
        }
        nw0 = nw1
      }
    }
  }

  def makect(nc: Int, c: Array[Double], startc: Int, ip: Array[Int]): Unit = {
    var j = 0
    var nch = 0
    var delta = .0
    var deltaj = .0
    ip(1) = nc
    if (nc > 1) {
      nch = nc >> 1
      delta = 0.785398163397448278999490867136046290 / nch
      c(startc) = cos(delta * nch)
      c(startc + nch) = 0.5 * c(startc)
      j = 1
      while ( {
        j < nch
      }) {
        deltaj = delta * j
        c(startc + j) = 0.5 * cos(deltaj)
        c(startc + nc - j) = 0.5 * sin(deltaj)

        j += 1
      }
    }
  }

  def makect(nc: Int, c: Array[Float], startc: Int, ip: Array[Int]): Unit = {
    var j = 0
    var nch = 0
    var delta = .0
    var deltaj = .0
    ip(1) = nc
    if (nc > 1) {
      nch = nc >> 1
      delta = 0.785398163397448278999490867136046290f / nch
      c(startc) = cos(delta * nch).toFloat
      c(startc + nch) = 0.5f * c(startc)
      j = 1
      while ( {
        j < nch
      }) {
        deltaj = delta * j
        c(startc + j) = 0.5f * cos(deltaj).toFloat
        c(startc + nc - j) = 0.5f * sin(deltaj).toFloat

        j += 1
      }
    }
  }

  def makewt(nw: Int, ip: Array[Int], w: Array[Float]): Unit = {
    var j       = 0
    var nwh     = 0
    var nw0     = 0
    var nw1     = 0
    var delta   = 0.0f
    var wn4r    = 0.0f
    var wk1r    = 0.0f
    var wk1i    = 0.0f
    var wk3r    = 0.0f
    var wk3i    = 0.0f
    var delta2  = 0.0f
    var deltaj  = 0.0f
    var deltaj3 = 0.0f
    ip(0) = nw
    ip(1) = 1
    if (nw > 2) {
      nwh = nw >> 1
      delta = 0.785398163397448278999490867136046290f / nwh
      delta2 = delta * 2
      wn4r = cos(delta * nwh).toFloat
      w(0) = 1
      w(1) = wn4r
      if (nwh == 4) {
        w(2) = cos(delta2).toFloat
        w(3) = sin(delta2).toFloat
      }
      else if (nwh > 4) {
        CommonUtils.makeipt(nw, ip)
        w(2) = 0.5f / cos(delta2).toFloat
        w(3) = 0.5f / cos(delta * 6).toFloat
        j = 4
        while ( {
          j < nwh
        }) {
          deltaj = delta * j
          deltaj3 = 3 * deltaj
          w(j) = cos(deltaj).toFloat
          w(j + 1) = sin(deltaj).toFloat
          w(j + 2) = cos(deltaj3).toFloat
          w(j + 3) = -sin(deltaj3).toFloat

          j += 4
        }
      }
      nw0 = 0
      while ( {
        nwh > 2
      }) {
        nw1 = nw0 + nwh
        nwh >>= 1
        w(nw1) = 1
        w(nw1 + 1) = wn4r
        if (nwh == 4) {
          wk1r = w(nw0 + 4)
          wk1i = w(nw0 + 5)
          w(nw1 + 2) = wk1r
          w(nw1 + 3) = wk1i
        }
        else if (nwh > 4) {
          wk1r = w(nw0 + 4)
          wk3r = w(nw0 + 6)
          w(nw1 + 2) = 0.5f / wk1r
          w(nw1 + 3) = 0.5f / wk3r
          j = 4
          while ( {
            j < nwh
          }) {
            val idx1 = nw0 + 2 * j
            val idx2 = nw1 + j
            wk1r = w(idx1)
            wk1i = w(idx1 + 1)
            wk3r = w(idx1 + 2)
            wk3i = w(idx1 + 3)
            w(idx2) = wk1r
            w(idx2 + 1) = wk1i
            w(idx2 + 2) = wk3r
            w(idx2 + 3) = wk3i

            j += 4
          }
        }
        nw0 = nw1
      }
    }
  }

  def cftfsub(n: Int, a: Array[Double], offa: Int, ip: Array[Int], nw: Int, w: Array[Double]): Unit = {
    if (n > 8) if (n > 32) {
      cftf1st(n, a, offa, w, nw - (n >> 2))
      if ((ConcurrencyUtils.numThreads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) cftrec4_th(n, a, offa, nw, w)
      else if (n > 512) cftrec4(n, a, offa, nw, w)
      else if (n > 128) cftleaf(n, 1, a, offa, nw, w)
      else cftfx41(n, a, offa, nw, w)
      bitrv2(n, ip, a, offa)
    }
    else if (n == 32) {
      cftf161(a, offa, w, nw - 8)
      bitrv216(a, offa)
    }
    else {
      cftf081(a, offa, w, 0)
      bitrv208(a, offa)
    }
    else if (n == 8) cftf040(a, offa)
    else if (n == 4) cftxb020(a, offa)
  }

  def cftbsub(n: Int, a: Array[Double], offa: Int, ip: Array[Int], nw: Int, w: Array[Double]): Unit = {
    if (n > 8) if (n > 32) {
      cftb1st(n, a, offa, w, nw - (n >> 2))
      if ((ConcurrencyUtils.numThreads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) cftrec4_th(n, a, offa, nw, w)
      else if (n > 512) cftrec4(n, a, offa, nw, w)
      else if (n > 128) cftleaf(n, 1, a, offa, nw, w)
      else cftfx41(n, a, offa, nw, w)
      bitrv2conj(n, ip, a, offa)
    }
    else if (n == 32) {
      cftf161(a, offa, w, nw - 8)
      bitrv216neg(a, offa)
    }
    else {
      cftf081(a, offa, w, 0)
      bitrv208neg(a, offa)
    }
    else if (n == 8) cftb040(a, offa)
    else if (n == 4) cftxb020(a, offa)
  }

  def bitrv2(n: Int, ip: Array[Int], a: Array[Double], offa: Int): Unit = {
    var j1 = 0
    var k1 = 0
    var l = 0
    var m = 0
    var nh = 0
    var nm = 0
    var xr = .0
    var xi = .0
    var yr = .0
    var yi = .0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    m = 1
    l = n >> 2
    while ( {
      l > 8
    }) {
      m <<= 1

      l >>= 2
    }
    nh = n >> 1
    nm = 4 * m
    if (l == 8) for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + 2 * ip(m + k)
        k1 = idx0 + 2 * ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + 2 * ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 += 2 * nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 -= nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= 2
      k1 -= nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nh + 2
      k1 += nh + 2
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= nh - nm
      k1 += 2 * nm - 2
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
    }
    else for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + ip(m + k)
        k1 = idx0 + ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 += nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
    }
  }

  def bitrv2conj(n: Int, ip: Array[Int], a: Array[Double], offa: Int): Unit = {
    var j1 = 0
    var k1 = 0
    var l = 0
    var m = 0
    var nh = 0
    var nm = 0
    var xr = .0
    var xi = .0
    var yr = .0
    var yi = .0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    m = 1
    l = n >> 2
    while ( {
      l > 8
    }) {
      m <<= 1

      l >>= 2
    }
    nh = n >> 1
    nm = 4 * m
    if (l == 8) for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + 2 * ip(m + k)
        k1 = idx0 + 2 * ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + 2 * ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
      j1 += nm
      k1 += 2 * nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 -= nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= 2
      k1 -= nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nh + 2
      k1 += nh + 2
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= nh - nm
      k1 += 2 * nm - 2
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
    }
    else for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + ip(m + k)
        k1 = idx0 + ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
      j1 += nm
      k1 += nm
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
    }
  }

  def bitrv216(a: Array[Double], offa: Int): Unit = {
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var x4r = .0
    var x4i = .0
    var x5r = .0
    var x5i = .0
    var x7r = .0
    var x7i = .0
    var x8r = .0
    var x8i = .0
    var x10r = .0
    var x10i = .0
    var x11r = .0
    var x11i = .0
    var x12r = .0
    var x12i = .0
    var x13r = .0
    var x13i = .0
    var x14r = .0
    var x14i = .0
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x2r = a(offa + 4)
    x2i = a(offa + 5)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x5r = a(offa + 10)
    x5i = a(offa + 11)
    x7r = a(offa + 14)
    x7i = a(offa + 15)
    x8r = a(offa + 16)
    x8i = a(offa + 17)
    x10r = a(offa + 20)
    x10i = a(offa + 21)
    x11r = a(offa + 22)
    x11i = a(offa + 23)
    x12r = a(offa + 24)
    x12i = a(offa + 25)
    x13r = a(offa + 26)
    x13i = a(offa + 27)
    x14r = a(offa + 28)
    x14i = a(offa + 29)
    a(offa + 2) = x8r
    a(offa + 3) = x8i
    a(offa + 4) = x4r
    a(offa + 5) = x4i
    a(offa + 6) = x12r
    a(offa + 7) = x12i
    a(offa + 8) = x2r
    a(offa + 9) = x2i
    a(offa + 10) = x10r
    a(offa + 11) = x10i
    a(offa + 14) = x14r
    a(offa + 15) = x14i
    a(offa + 16) = x1r
    a(offa + 17) = x1i
    a(offa + 20) = x5r
    a(offa + 21) = x5i
    a(offa + 22) = x13r
    a(offa + 23) = x13i
    a(offa + 24) = x3r
    a(offa + 25) = x3i
    a(offa + 26) = x11r
    a(offa + 27) = x11i
    a(offa + 28) = x7r
    a(offa + 29) = x7i
  }

  def bitrv216neg(a: Array[Double], offa: Int): Unit = {
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var x4r = .0
    var x4i = .0
    var x5r = .0
    var x5i = .0
    var x6r = .0
    var x6i = .0
    var x7r = .0
    var x7i = .0
    var x8r = .0
    var x8i = .0
    var x9r = .0
    var x9i = .0
    var x10r = .0
    var x10i = .0
    var x11r = .0
    var x11i = .0
    var x12r = .0
    var x12i = .0
    var x13r = .0
    var x13i = .0
    var x14r = .0
    var x14i = .0
    var x15r = .0
    var x15i = .0
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x2r = a(offa + 4)
    x2i = a(offa + 5)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x5r = a(offa + 10)
    x5i = a(offa + 11)
    x6r = a(offa + 12)
    x6i = a(offa + 13)
    x7r = a(offa + 14)
    x7i = a(offa + 15)
    x8r = a(offa + 16)
    x8i = a(offa + 17)
    x9r = a(offa + 18)
    x9i = a(offa + 19)
    x10r = a(offa + 20)
    x10i = a(offa + 21)
    x11r = a(offa + 22)
    x11i = a(offa + 23)
    x12r = a(offa + 24)
    x12i = a(offa + 25)
    x13r = a(offa + 26)
    x13i = a(offa + 27)
    x14r = a(offa + 28)
    x14i = a(offa + 29)
    x15r = a(offa + 30)
    x15i = a(offa + 31)
    a(offa + 2) = x15r
    a(offa + 3) = x15i
    a(offa + 4) = x7r
    a(offa + 5) = x7i
    a(offa + 6) = x11r
    a(offa + 7) = x11i
    a(offa + 8) = x3r
    a(offa + 9) = x3i
    a(offa + 10) = x13r
    a(offa + 11) = x13i
    a(offa + 12) = x5r
    a(offa + 13) = x5i
    a(offa + 14) = x9r
    a(offa + 15) = x9i
    a(offa + 16) = x1r
    a(offa + 17) = x1i
    a(offa + 18) = x14r
    a(offa + 19) = x14i
    a(offa + 20) = x6r
    a(offa + 21) = x6i
    a(offa + 22) = x10r
    a(offa + 23) = x10i
    a(offa + 24) = x2r
    a(offa + 25) = x2i
    a(offa + 26) = x12r
    a(offa + 27) = x12i
    a(offa + 28) = x4r
    a(offa + 29) = x4i
    a(offa + 30) = x8r
    a(offa + 31) = x8i
  }

  def bitrv208(a: Array[Double], offa: Int): Unit = {
    var x1r = .0
    var x1i = .0
    var x3r = .0
    var x3i = .0
    var x4r = .0
    var x4i = .0
    var x6r = .0
    var x6i = .0
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x6r = a(offa + 12)
    x6i = a(offa + 13)
    a(offa + 2) = x4r
    a(offa + 3) = x4i
    a(offa + 6) = x6r
    a(offa + 7) = x6i
    a(offa + 8) = x1r
    a(offa + 9) = x1i
    a(offa + 12) = x3r
    a(offa + 13) = x3i
  }

  def bitrv208neg(a: Array[Double], offa: Int): Unit = {
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var x4r = .0
    var x4i = .0
    var x5r = .0
    var x5i = .0
    var x6r = .0
    var x6i = .0
    var x7r = .0
    var x7i = .0
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x2r = a(offa + 4)
    x2i = a(offa + 5)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x5r = a(offa + 10)
    x5i = a(offa + 11)
    x6r = a(offa + 12)
    x6i = a(offa + 13)
    x7r = a(offa + 14)
    x7i = a(offa + 15)
    a(offa + 2) = x7r
    a(offa + 3) = x7i
    a(offa + 4) = x3r
    a(offa + 5) = x3i
    a(offa + 6) = x5r
    a(offa + 7) = x5i
    a(offa + 8) = x1r
    a(offa + 9) = x1i
    a(offa + 10) = x6r
    a(offa + 11) = x6i
    a(offa + 12) = x2r
    a(offa + 13) = x2i
    a(offa + 14) = x4r
    a(offa + 15) = x4i
  }

  def cftf1st(n: Int, a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var j0 = 0
    var j1 = 0
    var j2 = 0
    var j3 = 0
    var k = 0
    var m = 0
    var mh = 0
    var wn4r = .0
    var csc1 = .0
    var csc3 = .0
    var wk1r = .0
    var wk1i = .0
    var wk3r = .0
    var wk3i = .0
    var wd1r = .0
    var wd1i = .0
    var wd3r = .0
    var wd3i = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var y0r = .0
    var y0i = .0
    var y1r = .0
    var y1i = .0
    var y2r = .0
    var y2i = .0
    var y3r = .0
    var y3i = .0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    var idx3 = 0
    var idx4 = 0
    var idx5 = 0
    mh = n >> 3
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) + a(idx2)
    x0i = a(offa + 1) + a(idx2 + 1)
    x1r = a(offa) - a(idx2)
    x1i = a(offa + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    a(idx2) = x1r - x3i
    a(idx2 + 1) = x1i + x3r
    a(idx3) = x1r + x3i
    a(idx3 + 1) = x1i - x3r
    wn4r = w(startw + 1)
    csc1 = w(startw + 2)
    csc3 = w(startw + 3)
    wd1r = 1
    wd1i = 0
    wd3r = 1
    wd3i = 0
    k = 0
    var j = 2
    while ( {
      j < mh - 2
    }) {
      k += 4
      idx4 = startw + k
      wk1r = csc1 * (wd1r + w(idx4))
      wk1i = csc1 * (wd1i + w(idx4 + 1))
      wk3r = csc3 * (wd3r + w(idx4 + 2))
      wk3i = csc3 * (wd3i + w(idx4 + 3))
      wd1r = w(idx4)
      wd1i = w(idx4 + 1)
      wd3r = w(idx4 + 2)
      wd3i = w(idx4 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx5 = offa + j
      x0r = a(idx5) + a(idx2)
      x0i = a(idx5 + 1) + a(idx2 + 1)
      x1r = a(idx5) - a(idx2)
      x1i = a(idx5 + 1) - a(idx2 + 1)
      y0r = a(idx5 + 2) + a(idx2 + 2)
      y0i = a(idx5 + 3) + a(idx2 + 3)
      y1r = a(idx5 + 2) - a(idx2 + 2)
      y1i = a(idx5 + 3) - a(idx2 + 3)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 + 2) + a(idx3 + 2)
      y2i = a(idx1 + 3) + a(idx3 + 3)
      y3r = a(idx1 + 2) - a(idx3 + 2)
      y3i = a(idx1 + 3) - a(idx3 + 3)
      a(idx5) = x0r + x2r
      a(idx5 + 1) = x0i + x2i
      a(idx5 + 2) = y0r + y2r
      a(idx5 + 3) = y0i + y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      a(idx1 + 2) = y0r - y2r
      a(idx1 + 3) = y0i - y2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1r * x0r - wk1i * x0i
      a(idx2 + 1) = wk1r * x0i + wk1i * x0r
      x0r = y1r - y3i
      x0i = y1i + y3r
      a(idx2 + 2) = wd1r * x0r - wd1i * x0i
      a(idx2 + 3) = wd1r * x0i + wd1i * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3r * x0r + wk3i * x0i
      a(idx3 + 1) = wk3r * x0i - wk3i * x0r
      x0r = y1r + y3i
      x0i = y1i - y3r
      a(idx3 + 2) = wd3r * x0r + wd3i * x0i
      a(idx3 + 3) = wd3r * x0i - wd3i * x0r
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) + a(idx2)
      x0i = a(idx0 + 1) + a(idx2 + 1)
      x1r = a(idx0) - a(idx2)
      x1i = a(idx0 + 1) - a(idx2 + 1)
      y0r = a(idx0 - 2) + a(idx2 - 2)
      y0i = a(idx0 - 1) + a(idx2 - 1)
      y1r = a(idx0 - 2) - a(idx2 - 2)
      y1i = a(idx0 - 1) - a(idx2 - 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 - 2) + a(idx3 - 2)
      y2i = a(idx1 - 1) + a(idx3 - 1)
      y3r = a(idx1 - 2) - a(idx3 - 2)
      y3i = a(idx1 - 1) - a(idx3 - 1)
      a(idx0) = x0r + x2r
      a(idx0 + 1) = x0i + x2i
      a(idx0 - 2) = y0r + y2r
      a(idx0 - 1) = y0i + y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      a(idx1 - 2) = y0r - y2r
      a(idx1 - 1) = y0i - y2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1i * x0r - wk1r * x0i
      a(idx2 + 1) = wk1i * x0i + wk1r * x0r
      x0r = y1r - y3i
      x0i = y1i + y3r
      a(idx2 - 2) = wd1i * x0r - wd1r * x0i
      a(idx2 - 1) = wd1i * x0i + wd1r * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3i * x0r + wk3r * x0i
      a(idx3 + 1) = wk3i * x0i - wk3r * x0r
      x0r = y1r + y3i
      x0i = y1i - y3r
      a(offa + j3 - 2) = wd3i * x0r + wd3r * x0i
      a(offa + j3 - 1) = wd3i * x0i - wd3r * x0r

      j += 4
    }
    wk1r = csc1 * (wd1r + wn4r)
    wk1i = csc1 * (wd1i + wn4r)
    wk3r = csc3 * (wd3r - wn4r)
    wk3i = csc3 * (wd3i - wn4r)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0 - 2) + a(idx2 - 2)
    x0i = a(idx0 - 1) + a(idx2 - 1)
    x1r = a(idx0 - 2) - a(idx2 - 2)
    x1i = a(idx0 - 1) - a(idx2 - 1)
    x2r = a(idx1 - 2) + a(idx3 - 2)
    x2i = a(idx1 - 1) + a(idx3 - 1)
    x3r = a(idx1 - 2) - a(idx3 - 2)
    x3i = a(idx1 - 1) - a(idx3 - 1)
    a(idx0 - 2) = x0r + x2r
    a(idx0 - 1) = x0i + x2i
    a(idx1 - 2) = x0r - x2r
    a(idx1 - 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2 - 2) = wk1r * x0r - wk1i * x0i
    a(idx2 - 1) = wk1r * x0i + wk1i * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3 - 2) = wk3r * x0r + wk3i * x0i
    a(idx3 - 1) = wk3r * x0i - wk3i * x0r
    x0r = a(idx0) + a(idx2)
    x0i = a(idx0 + 1) + a(idx2 + 1)
    x1r = a(idx0) - a(idx2)
    x1i = a(idx0 + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(idx0) = x0r + x2r
    a(idx0 + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2) = wn4r * (x0r - x0i)
    a(idx2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3) = -wn4r * (x0r + x0i)
    a(idx3 + 1) = -wn4r * (x0i - x0r)
    x0r = a(idx0 + 2) + a(idx2 + 2)
    x0i = a(idx0 + 3) + a(idx2 + 3)
    x1r = a(idx0 + 2) - a(idx2 + 2)
    x1i = a(idx0 + 3) - a(idx2 + 3)
    x2r = a(idx1 + 2) + a(idx3 + 2)
    x2i = a(idx1 + 3) + a(idx3 + 3)
    x3r = a(idx1 + 2) - a(idx3 + 2)
    x3i = a(idx1 + 3) - a(idx3 + 3)
    a(idx0 + 2) = x0r + x2r
    a(idx0 + 3) = x0i + x2i
    a(idx1 + 2) = x0r - x2r
    a(idx1 + 3) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2 + 2) = wk1i * x0r - wk1r * x0i
    a(idx2 + 3) = wk1i * x0i + wk1r * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3 + 2) = wk3i * x0r + wk3r * x0i
    a(idx3 + 3) = wk3i * x0i - wk3r * x0r
  }

  def cftb1st(n: Int, a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var j0 = 0
    var j1 = 0
    var j2 = 0
    var j3 = 0
    var k = 0
    var m = 0
    var mh = 0
    var wn4r = .0
    var csc1 = .0
    var csc3 = .0
    var wk1r = .0
    var wk1i = .0
    var wk3r = .0
    var wk3i = .0
    var wd1r = .0
    var wd1i = .0
    var wd3r = .0
    var wd3i = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var y0r = .0
    var y0i = .0
    var y1r = .0
    var y1i = .0
    var y2r = .0
    var y2i = .0
    var y3r = .0
    var y3i = .0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    var idx3 = 0
    var idx4 = 0
    var idx5 = 0
    mh = n >> 3
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) + a(idx2)
    x0i = -a(offa + 1) - a(idx2 + 1)
    x1r = a(offa) - a(idx2)
    x1i = -a(offa + 1) + a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i - x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i + x2i
    a(idx2) = x1r + x3i
    a(idx2 + 1) = x1i + x3r
    a(idx3) = x1r - x3i
    a(idx3 + 1) = x1i - x3r
    wn4r = w(startw + 1)
    csc1 = w(startw + 2)
    csc3 = w(startw + 3)
    wd1r = 1
    wd1i = 0
    wd3r = 1
    wd3i = 0
    k = 0
    var j = 2
    while ( {
      j < mh - 2
    }) {
      k += 4
      idx4 = startw + k
      wk1r = csc1 * (wd1r + w(idx4))
      wk1i = csc1 * (wd1i + w(idx4 + 1))
      wk3r = csc3 * (wd3r + w(idx4 + 2))
      wk3i = csc3 * (wd3i + w(idx4 + 3))
      wd1r = w(idx4)
      wd1i = w(idx4 + 1)
      wd3r = w(idx4 + 2)
      wd3i = w(idx4 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx5 = offa + j
      x0r = a(idx5) + a(idx2)
      x0i = -a(idx5 + 1) - a(idx2 + 1)
      x1r = a(idx5) - a(offa + j2)
      x1i = -a(idx5 + 1) + a(idx2 + 1)
      y0r = a(idx5 + 2) + a(idx2 + 2)
      y0i = -a(idx5 + 3) - a(idx2 + 3)
      y1r = a(idx5 + 2) - a(idx2 + 2)
      y1i = -a(idx5 + 3) + a(idx2 + 3)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 + 2) + a(idx3 + 2)
      y2i = a(idx1 + 3) + a(idx3 + 3)
      y3r = a(idx1 + 2) - a(idx3 + 2)
      y3i = a(idx1 + 3) - a(idx3 + 3)
      a(idx5) = x0r + x2r
      a(idx5 + 1) = x0i - x2i
      a(idx5 + 2) = y0r + y2r
      a(idx5 + 3) = y0i - y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i + x2i
      a(idx1 + 2) = y0r - y2r
      a(idx1 + 3) = y0i + y2i
      x0r = x1r + x3i
      x0i = x1i + x3r
      a(idx2) = wk1r * x0r - wk1i * x0i
      a(idx2 + 1) = wk1r * x0i + wk1i * x0r
      x0r = y1r + y3i
      x0i = y1i + y3r
      a(idx2 + 2) = wd1r * x0r - wd1i * x0i
      a(idx2 + 3) = wd1r * x0i + wd1i * x0r
      x0r = x1r - x3i
      x0i = x1i - x3r
      a(idx3) = wk3r * x0r + wk3i * x0i
      a(idx3 + 1) = wk3r * x0i - wk3i * x0r
      x0r = y1r - y3i
      x0i = y1i - y3r
      a(idx3 + 2) = wd3r * x0r + wd3i * x0i
      a(idx3 + 3) = wd3r * x0i - wd3i * x0r
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) + a(idx2)
      x0i = -a(idx0 + 1) - a(idx2 + 1)
      x1r = a(idx0) - a(idx2)
      x1i = -a(idx0 + 1) + a(idx2 + 1)
      y0r = a(idx0 - 2) + a(idx2 - 2)
      y0i = -a(idx0 - 1) - a(idx2 - 1)
      y1r = a(idx0 - 2) - a(idx2 - 2)
      y1i = -a(idx0 - 1) + a(idx2 - 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 - 2) + a(idx3 - 2)
      y2i = a(idx1 - 1) + a(idx3 - 1)
      y3r = a(idx1 - 2) - a(idx3 - 2)
      y3i = a(idx1 - 1) - a(idx3 - 1)
      a(idx0) = x0r + x2r
      a(idx0 + 1) = x0i - x2i
      a(idx0 - 2) = y0r + y2r
      a(idx0 - 1) = y0i - y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i + x2i
      a(idx1 - 2) = y0r - y2r
      a(idx1 - 1) = y0i + y2i
      x0r = x1r + x3i
      x0i = x1i + x3r
      a(idx2) = wk1i * x0r - wk1r * x0i
      a(idx2 + 1) = wk1i * x0i + wk1r * x0r
      x0r = y1r + y3i
      x0i = y1i + y3r
      a(idx2 - 2) = wd1i * x0r - wd1r * x0i
      a(idx2 - 1) = wd1i * x0i + wd1r * x0r
      x0r = x1r - x3i
      x0i = x1i - x3r
      a(idx3) = wk3i * x0r + wk3r * x0i
      a(idx3 + 1) = wk3i * x0i - wk3r * x0r
      x0r = y1r - y3i
      x0i = y1i - y3r
      a(idx3 - 2) = wd3i * x0r + wd3r * x0i
      a(idx3 - 1) = wd3i * x0i - wd3r * x0r

      j += 4
    }
    wk1r = csc1 * (wd1r + wn4r)
    wk1i = csc1 * (wd1i + wn4r)
    wk3r = csc3 * (wd3r - wn4r)
    wk3i = csc3 * (wd3i - wn4r)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0 - 2) + a(idx2 - 2)
    x0i = -a(idx0 - 1) - a(idx2 - 1)
    x1r = a(idx0 - 2) - a(idx2 - 2)
    x1i = -a(idx0 - 1) + a(idx2 - 1)
    x2r = a(idx1 - 2) + a(idx3 - 2)
    x2i = a(idx1 - 1) + a(idx3 - 1)
    x3r = a(idx1 - 2) - a(idx3 - 2)
    x3i = a(idx1 - 1) - a(idx3 - 1)
    a(idx0 - 2) = x0r + x2r
    a(idx0 - 1) = x0i - x2i
    a(idx1 - 2) = x0r - x2r
    a(idx1 - 1) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(idx2 - 2) = wk1r * x0r - wk1i * x0i
    a(idx2 - 1) = wk1r * x0i + wk1i * x0r
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(idx3 - 2) = wk3r * x0r + wk3i * x0i
    a(idx3 - 1) = wk3r * x0i - wk3i * x0r
    x0r = a(idx0) + a(idx2)
    x0i = -a(idx0 + 1) - a(idx2 + 1)
    x1r = a(idx0) - a(idx2)
    x1i = -a(idx0 + 1) + a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(idx0) = x0r + x2r
    a(idx0 + 1) = x0i - x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(idx2) = wn4r * (x0r - x0i)
    a(idx2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(idx3) = -wn4r * (x0r + x0i)
    a(idx3 + 1) = -wn4r * (x0i - x0r)
    x0r = a(idx0 + 2) + a(idx2 + 2)
    x0i = -a(idx0 + 3) - a(idx2 + 3)
    x1r = a(idx0 + 2) - a(idx2 + 2)
    x1i = -a(idx0 + 3) + a(idx2 + 3)
    x2r = a(idx1 + 2) + a(idx3 + 2)
    x2i = a(idx1 + 3) + a(idx3 + 3)
    x3r = a(idx1 + 2) - a(idx3 + 2)
    x3i = a(idx1 + 3) - a(idx3 + 3)
    a(idx0 + 2) = x0r + x2r
    a(idx0 + 3) = x0i - x2i
    a(idx1 + 2) = x0r - x2r
    a(idx1 + 3) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(idx2 + 2) = wk1i * x0r - wk1r * x0i
    a(idx2 + 3) = wk1i * x0i + wk1r * x0r
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(idx3 + 2) = wk3i * x0r + wk3r * x0i
    a(idx3 + 3) = wk3i * x0i - wk3r * x0r
  }

  def cftrec4_th(n: Int, a: Array[Double], offa: Int, nw: Int, w: Array[Double]): Unit = {
    var i         = 0
    var idiv4     = 0
    var m         = 0
    var nthreads  = 0
    var idx       = 0
    nthreads      = 2
    idiv4         = 0
    m = n >> 1
    if (n >= CommonUtils.threadsBeginN_1D_FFT_4Threads) {
      nthreads = 4
      idiv4 = 1
      m >>= 1
    }
    val futures = new Array[Future[_]](nthreads)
    val mf = m
    i = 0
    while (i < nthreads) {
      val firstIdx = offa + i * m
      futures(idx) = if (i != idiv4) Future {
        var isplt = 0
        var j     = 0
        var k     = 0
        var m     = 0
        val idx1  = firstIdx + mf
        m = n
        while (m > 512) {
          m >>= 2
          cftmdl1(m, a, idx1 - m, w, nw - (m >> 1))
        }
        cftleaf(m, 1, a, idx1 - m, nw, w)
        k = 0
        val idx2 = firstIdx - m
        j = mf - m
        while (j > 0) {
          k += 1
          isplt = cfttree(m, j, k, a, firstIdx, nw, w)
          cftleaf(m, isplt, a, idx2 + j, nw, w)

          j -= m
        }
      } else Future {
        var isplt = 0
        var j     = 0
        var k     = 0
        var m     = 0
        val idx1  = firstIdx + mf
        k = 1
        m = n
        while (m > 512) {
          m >>= 2
          k <<= 2
          cftmdl2(m, a, idx1 - m, w, nw - m)
        }
        cftleaf(m, 0, a, idx1 - m, nw, w)
        k >>= 1
        val idx2 = firstIdx - m
        j = mf - m
        while (j > 0) {
          k += 1
          isplt = cfttree(m, j, k, a, firstIdx, nw, w)
          cftleaf(m, isplt, a, idx2 + j, nw, w)

          j -= m
        }
      }
      idx += 1

      i += 1
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  private final val name = "CommonUtils"

  def cftrec4(n: Int, a: Array[Double], offa: Int, nw: Int, w: Array[Double]): Unit = {
    var isplt = 0
    var j = 0
    var k = 0
    var m = 0
    m = n
    val idx1 = offa + n
    while ( {
      m > 512
    }) {
      m >>= 2
      cftmdl1(m, a, idx1 - m, w, nw - (m >> 1))
    }
    cftleaf(m, 1, a, idx1 - m, nw, w)
    k = 0
    val idx2 = offa - m
    j = n - m
    while ( {
      j > 0
    }) {
      k += 1
      isplt = cfttree(m, j, k, a, offa, nw, w)
      cftleaf(m, isplt, a, idx2 + j, nw, w)

      j -= m
    }
  }

  def cfttree(n: Int, j: Int, k: Int, a: Array[Double], offa: Int, nw: Int, w: Array[Double]): Int = {
    var i = 0
    var isplt = 0
    var m = 0
    val idx1 = offa - n
    if ((k & 3) != 0) {
      isplt = k & 1
      if (isplt != 0) cftmdl1(n, a, idx1 + j, w, nw - (n >> 1))
      else cftmdl2(n, a, idx1 + j, w, nw - n)
    }
    else {
      m = n
      i = k
      while ( {
        (i & 3) == 0
      }) {
        m <<= 2

        i >>= 2
      }
      isplt = i & 1
      val idx2 = offa + j
      if (isplt != 0) while ( {
        m > 128
      }) {
        cftmdl1(m, a, idx2 - m, w, nw - (m >> 1))
        m >>= 2
      }
      else while ( {
        m > 128
      }) {
        cftmdl2(m, a, idx2 - m, w, nw - m)
        m >>= 2
      }
    }
    isplt
  }

  def cftleaf(n: Int, isplt: Int, a: Array[Double], offa: Int, nw: Int, w: Array[Double]): Unit = {
    if (n == 512) {
      cftmdl1(128, a, offa, w, nw - 64)
      cftf161(a, offa, w, nw - 8)
      cftf162(a, offa + 32, w, nw - 32)
      cftf161(a, offa + 64, w, nw - 8)
      cftf161(a, offa + 96, w, nw - 8)
      cftmdl2(128, a, offa + 128, w, nw - 128)
      cftf161(a, offa + 128, w, nw - 8)
      cftf162(a, offa + 160, w, nw - 32)
      cftf161(a, offa + 192, w, nw - 8)
      cftf162(a, offa + 224, w, nw - 32)
      cftmdl1(128, a, offa + 256, w, nw - 64)
      cftf161(a, offa + 256, w, nw - 8)
      cftf162(a, offa + 288, w, nw - 32)
      cftf161(a, offa + 320, w, nw - 8)
      cftf161(a, offa + 352, w, nw - 8)
      if (isplt != 0) {
        cftmdl1(128, a, offa + 384, w, nw - 64)
        cftf161(a, offa + 480, w, nw - 8)
      }
      else {
        cftmdl2(128, a, offa + 384, w, nw - 128)
        cftf162(a, offa + 480, w, nw - 32)
      }
      cftf161(a, offa + 384, w, nw - 8)
      cftf162(a, offa + 416, w, nw - 32)
      cftf161(a, offa + 448, w, nw - 8)
    }
    else {
      cftmdl1(64, a, offa, w, nw - 32)
      cftf081(a, offa, w, nw - 8)
      cftf082(a, offa + 16, w, nw - 8)
      cftf081(a, offa + 32, w, nw - 8)
      cftf081(a, offa + 48, w, nw - 8)
      cftmdl2(64, a, offa + 64, w, nw - 64)
      cftf081(a, offa + 64, w, nw - 8)
      cftf082(a, offa + 80, w, nw - 8)
      cftf081(a, offa + 96, w, nw - 8)
      cftf082(a, offa + 112, w, nw - 8)
      cftmdl1(64, a, offa + 128, w, nw - 32)
      cftf081(a, offa + 128, w, nw - 8)
      cftf082(a, offa + 144, w, nw - 8)
      cftf081(a, offa + 160, w, nw - 8)
      cftf081(a, offa + 176, w, nw - 8)
      if (isplt != 0) {
        cftmdl1(64, a, offa + 192, w, nw - 32)
        cftf081(a, offa + 240, w, nw - 8)
      }
      else {
        cftmdl2(64, a, offa + 192, w, nw - 64)
        cftf082(a, offa + 240, w, nw - 8)
      }
      cftf081(a, offa + 192, w, nw - 8)
      cftf082(a, offa + 208, w, nw - 8)
      cftf081(a, offa + 224, w, nw - 8)
    }
  }

  def cftmdl1(n: Int, a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var j0 = 0
    var j1 = 0
    var j2 = 0
    var j3 = 0
    var k = 0
    var m = 0
    var mh = 0
    var wn4r = .0
    var wk1r = .0
    var wk1i = .0
    var wk3r = .0
    var wk3i = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    var idx3 = 0
    var idx4 = 0
    var idx5 = 0
    mh = n >> 3
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) + a(idx2)
    x0i = a(offa + 1) + a(idx2 + 1)
    x1r = a(offa) - a(idx2)
    x1i = a(offa + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    a(idx2) = x1r - x3i
    a(idx2 + 1) = x1i + x3r
    a(idx3) = x1r + x3i
    a(idx3 + 1) = x1i - x3r
    wn4r = w(startw + 1)
    k = 0
    var j = 2
    while ( {
      j < mh
    }) {
      k += 4
      idx4 = startw + k
      wk1r = w(idx4)
      wk1i = w(idx4 + 1)
      wk3r = w(idx4 + 2)
      wk3i = w(idx4 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx5 = offa + j
      x0r = a(idx5) + a(idx2)
      x0i = a(idx5 + 1) + a(idx2 + 1)
      x1r = a(idx5) - a(idx2)
      x1i = a(idx5 + 1) - a(idx2 + 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      a(idx5) = x0r + x2r
      a(idx5 + 1) = x0i + x2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1r * x0r - wk1i * x0i
      a(idx2 + 1) = wk1r * x0i + wk1i * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3r * x0r + wk3i * x0i
      a(idx3 + 1) = wk3r * x0i - wk3i * x0r
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) + a(idx2)
      x0i = a(idx0 + 1) + a(idx2 + 1)
      x1r = a(idx0) - a(idx2)
      x1i = a(idx0 + 1) - a(idx2 + 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      a(idx0) = x0r + x2r
      a(idx0 + 1) = x0i + x2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1i * x0r - wk1r * x0i
      a(idx2 + 1) = wk1i * x0i + wk1r * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3i * x0r + wk3r * x0i
      a(idx3 + 1) = wk3i * x0i - wk3r * x0r

      j += 2
    }
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0) + a(idx2)
    x0i = a(idx0 + 1) + a(idx2 + 1)
    x1r = a(idx0) - a(idx2)
    x1i = a(idx0 + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(idx0) = x0r + x2r
    a(idx0 + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2) = wn4r * (x0r - x0i)
    a(idx2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3) = -wn4r * (x0r + x0i)
    a(idx3 + 1) = -wn4r * (x0i - x0r)
  }

  def cftmdl2(n: Int, a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var j0 = 0
    var j1 = 0
    var j2 = 0
    var j3 = 0
    var k = 0
    var kr = 0
    var m = 0
    var mh = 0
    var wn4r = .0
    var wk1r = .0
    var wk1i = .0
    var wk3r = .0
    var wk3i = .0
    var wd1r = .0
    var wd1i = .0
    var wd3r = .0
    var wd3i = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var y0r = .0
    var y0i = .0
    var y2r = .0
    var y2i = .0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    var idx3 = 0
    var idx4 = 0
    var idx5 = 0
    var idx6 = 0
    mh = n >> 3
    m = 2 * mh
    wn4r = w(startw + 1)
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) - a(idx2 + 1)
    x0i = a(offa + 1) + a(idx2)
    x1r = a(offa) + a(idx2 + 1)
    x1i = a(offa + 1) - a(idx2)
    x2r = a(idx1) - a(idx3 + 1)
    x2i = a(idx1 + 1) + a(idx3)
    x3r = a(idx1) + a(idx3 + 1)
    x3i = a(idx1 + 1) - a(idx3)
    y0r = wn4r * (x2r - x2i)
    y0i = wn4r * (x2i + x2r)
    a(offa) = x0r + y0r
    a(offa + 1) = x0i + y0i
    a(idx1) = x0r - y0r
    a(idx1 + 1) = x0i - y0i
    y0r = wn4r * (x3r - x3i)
    y0i = wn4r * (x3i + x3r)
    a(idx2) = x1r - y0i
    a(idx2 + 1) = x1i + y0r
    a(idx3) = x1r + y0i
    a(idx3 + 1) = x1i - y0r
    k = 0
    kr = 2 * m
    var j = 2
    while ( {
      j < mh
    }) {
      k += 4
      idx4 = startw + k
      wk1r = w(idx4)
      wk1i = w(idx4 + 1)
      wk3r = w(idx4 + 2)
      wk3i = w(idx4 + 3)
      kr -= 4
      idx5 = startw + kr
      wd1i = w(idx5)
      wd1r = w(idx5 + 1)
      wd3i = w(idx5 + 2)
      wd3r = w(idx5 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx6 = offa + j
      x0r = a(idx6) - a(idx2 + 1)
      x0i = a(idx6 + 1) + a(idx2)
      x1r = a(idx6) + a(idx2 + 1)
      x1i = a(idx6 + 1) - a(idx2)
      x2r = a(idx1) - a(idx3 + 1)
      x2i = a(idx1 + 1) + a(idx3)
      x3r = a(idx1) + a(idx3 + 1)
      x3i = a(idx1 + 1) - a(idx3)
      y0r = wk1r * x0r - wk1i * x0i
      y0i = wk1r * x0i + wk1i * x0r
      y2r = wd1r * x2r - wd1i * x2i
      y2i = wd1r * x2i + wd1i * x2r
      a(idx6) = y0r + y2r
      a(idx6 + 1) = y0i + y2i
      a(idx1) = y0r - y2r
      a(idx1 + 1) = y0i - y2i
      y0r = wk3r * x1r + wk3i * x1i
      y0i = wk3r * x1i - wk3i * x1r
      y2r = wd3r * x3r + wd3i * x3i
      y2i = wd3r * x3i - wd3i * x3r
      a(idx2) = y0r + y2r
      a(idx2 + 1) = y0i + y2i
      a(idx3) = y0r - y2r
      a(idx3 + 1) = y0i - y2i
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) - a(idx2 + 1)
      x0i = a(idx0 + 1) + a(idx2)
      x1r = a(idx0) + a(idx2 + 1)
      x1i = a(idx0 + 1) - a(idx2)
      x2r = a(idx1) - a(idx3 + 1)
      x2i = a(idx1 + 1) + a(idx3)
      x3r = a(idx1) + a(idx3 + 1)
      x3i = a(idx1 + 1) - a(idx3)
      y0r = wd1i * x0r - wd1r * x0i
      y0i = wd1i * x0i + wd1r * x0r
      y2r = wk1i * x2r - wk1r * x2i
      y2i = wk1i * x2i + wk1r * x2r
      a(idx0) = y0r + y2r
      a(idx0 + 1) = y0i + y2i
      a(idx1) = y0r - y2r
      a(idx1 + 1) = y0i - y2i
      y0r = wd3i * x1r + wd3r * x1i
      y0i = wd3i * x1i - wd3r * x1r
      y2r = wk3i * x3r + wk3r * x3i
      y2i = wk3i * x3i - wk3r * x3r
      a(idx2) = y0r + y2r
      a(idx2 + 1) = y0i + y2i
      a(idx3) = y0r - y2r
      a(idx3 + 1) = y0i - y2i

      j += 2
    }
    wk1r = w(startw + m)
    wk1i = w(startw + m + 1)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0) - a(idx2 + 1)
    x0i = a(idx0 + 1) + a(idx2)
    x1r = a(idx0) + a(idx2 + 1)
    x1i = a(idx0 + 1) - a(idx2)
    x2r = a(idx1) - a(idx3 + 1)
    x2i = a(idx1 + 1) + a(idx3)
    x3r = a(idx1) + a(idx3 + 1)
    x3i = a(idx1 + 1) - a(idx3)
    y0r = wk1r * x0r - wk1i * x0i
    y0i = wk1r * x0i + wk1i * x0r
    y2r = wk1i * x2r - wk1r * x2i
    y2i = wk1i * x2i + wk1r * x2r
    a(idx0) = y0r + y2r
    a(idx0 + 1) = y0i + y2i
    a(idx1) = y0r - y2r
    a(idx1 + 1) = y0i - y2i
    y0r = wk1i * x1r - wk1r * x1i
    y0i = wk1i * x1i + wk1r * x1r
    y2r = wk1r * x3r - wk1i * x3i
    y2i = wk1r * x3i + wk1i * x3r
    a(idx2) = y0r - y2r
    a(idx2 + 1) = y0i - y2i
    a(idx3) = y0r + y2r
    a(idx3 + 1) = y0i + y2i
  }

  def cftfx41(n: Int, a: Array[Double], offa: Int, nw: Int, w: Array[Double]): Unit = {
    if (n == 128) {
      cftf161(a, offa, w, nw - 8)
      cftf162(a, offa + 32, w, nw - 32)
      cftf161(a, offa + 64, w, nw - 8)
      cftf161(a, offa + 96, w, nw - 8)
    }
    else {
      cftf081(a, offa, w, nw - 8)
      cftf082(a, offa + 16, w, nw - 8)
      cftf081(a, offa + 32, w, nw - 8)
      cftf081(a, offa + 48, w, nw - 8)
    }
  }

  def cftf161(a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var wn4r = .0
    var wk1r = .0
    var wk1i = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var y0r = .0
    var y0i = .0
    var y1r = .0
    var y1i = .0
    var y2r = .0
    var y2i = .0
    var y3r = .0
    var y3i = .0
    var y4r = .0
    var y4i = .0
    var y5r = .0
    var y5i = .0
    var y6r = .0
    var y6i = .0
    var y7r = .0
    var y7i = .0
    var y8r = .0
    var y8i = .0
    var y9r = .0
    var y9i = .0
    var y10r = .0
    var y10i = .0
    var y11r = .0
    var y11i = .0
    var y12r = .0
    var y12i = .0
    var y13r = .0
    var y13i = .0
    var y14r = .0
    var y14i = .0
    var y15r = .0
    var y15i = .0
    wn4r = w(startw + 1)
    wk1r = w(startw + 2)
    wk1i = w(startw + 3)
    x0r = a(offa) + a(offa + 16)
    x0i = a(offa + 1) + a(offa + 17)
    x1r = a(offa) - a(offa + 16)
    x1i = a(offa + 1) - a(offa + 17)
    x2r = a(offa + 8) + a(offa + 24)
    x2i = a(offa + 9) + a(offa + 25)
    x3r = a(offa + 8) - a(offa + 24)
    x3i = a(offa + 9) - a(offa + 25)
    y0r = x0r + x2r
    y0i = x0i + x2i
    y4r = x0r - x2r
    y4i = x0i - x2i
    y8r = x1r - x3i
    y8i = x1i + x3r
    y12r = x1r + x3i
    y12i = x1i - x3r
    x0r = a(offa + 2) + a(offa + 18)
    x0i = a(offa + 3) + a(offa + 19)
    x1r = a(offa + 2) - a(offa + 18)
    x1i = a(offa + 3) - a(offa + 19)
    x2r = a(offa + 10) + a(offa + 26)
    x2i = a(offa + 11) + a(offa + 27)
    x3r = a(offa + 10) - a(offa + 26)
    x3i = a(offa + 11) - a(offa + 27)
    y1r = x0r + x2r
    y1i = x0i + x2i
    y5r = x0r - x2r
    y5i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y9r = wk1r * x0r - wk1i * x0i
    y9i = wk1r * x0i + wk1i * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    y13r = wk1i * x0r - wk1r * x0i
    y13i = wk1i * x0i + wk1r * x0r
    x0r = a(offa + 4) + a(offa + 20)
    x0i = a(offa + 5) + a(offa + 21)
    x1r = a(offa + 4) - a(offa + 20)
    x1i = a(offa + 5) - a(offa + 21)
    x2r = a(offa + 12) + a(offa + 28)
    x2i = a(offa + 13) + a(offa + 29)
    x3r = a(offa + 12) - a(offa + 28)
    x3i = a(offa + 13) - a(offa + 29)
    y2r = x0r + x2r
    y2i = x0i + x2i
    y6r = x0r - x2r
    y6i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y10r = wn4r * (x0r - x0i)
    y10i = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    y14r = wn4r * (x0r + x0i)
    y14i = wn4r * (x0i - x0r)
    x0r = a(offa + 6) + a(offa + 22)
    x0i = a(offa + 7) + a(offa + 23)
    x1r = a(offa + 6) - a(offa + 22)
    x1i = a(offa + 7) - a(offa + 23)
    x2r = a(offa + 14) + a(offa + 30)
    x2i = a(offa + 15) + a(offa + 31)
    x3r = a(offa + 14) - a(offa + 30)
    x3i = a(offa + 15) - a(offa + 31)
    y3r = x0r + x2r
    y3i = x0i + x2i
    y7r = x0r - x2r
    y7i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y11r = wk1i * x0r - wk1r * x0i
    y11i = wk1i * x0i + wk1r * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    y15r = wk1r * x0r - wk1i * x0i
    y15i = wk1r * x0i + wk1i * x0r
    x0r = y12r - y14r
    x0i = y12i - y14i
    x1r = y12r + y14r
    x1i = y12i + y14i
    x2r = y13r - y15r
    x2i = y13i - y15i
    x3r = y13r + y15r
    x3i = y13i + y15i
    a(offa + 24) = x0r + x2r
    a(offa + 25) = x0i + x2i
    a(offa + 26) = x0r - x2r
    a(offa + 27) = x0i - x2i
    a(offa + 28) = x1r - x3i
    a(offa + 29) = x1i + x3r
    a(offa + 30) = x1r + x3i
    a(offa + 31) = x1i - x3r
    x0r = y8r + y10r
    x0i = y8i + y10i
    x1r = y8r - y10r
    x1i = y8i - y10i
    x2r = y9r + y11r
    x2i = y9i + y11i
    x3r = y9r - y11r
    x3i = y9i - y11i
    a(offa + 16) = x0r + x2r
    a(offa + 17) = x0i + x2i
    a(offa + 18) = x0r - x2r
    a(offa + 19) = x0i - x2i
    a(offa + 20) = x1r - x3i
    a(offa + 21) = x1i + x3r
    a(offa + 22) = x1r + x3i
    a(offa + 23) = x1i - x3r
    x0r = y5r - y7i
    x0i = y5i + y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    x0r = y5r + y7i
    x0i = y5i - y7r
    x3r = wn4r * (x0r - x0i)
    x3i = wn4r * (x0i + x0r)
    x0r = y4r - y6i
    x0i = y4i + y6r
    x1r = y4r + y6i
    x1i = y4i - y6r
    a(offa + 8) = x0r + x2r
    a(offa + 9) = x0i + x2i
    a(offa + 10) = x0r - x2r
    a(offa + 11) = x0i - x2i
    a(offa + 12) = x1r - x3i
    a(offa + 13) = x1i + x3r
    a(offa + 14) = x1r + x3i
    a(offa + 15) = x1i - x3r
    x0r = y0r + y2r
    x0i = y0i + y2i
    x1r = y0r - y2r
    x1i = y0i - y2i
    x2r = y1r + y3r
    x2i = y1i + y3i
    x3r = y1r - y3r
    x3i = y1i - y3i
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(offa + 2) = x0r - x2r
    a(offa + 3) = x0i - x2i
    a(offa + 4) = x1r - x3i
    a(offa + 5) = x1i + x3r
    a(offa + 6) = x1r + x3i
    a(offa + 7) = x1i - x3r
  }

  def cftf162(a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var wn4r = .0
    var wk1r = .0
    var wk1i = .0
    var wk2r = .0
    var wk2i = .0
    var wk3r = .0
    var wk3i = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var y0r = .0
    var y0i = .0
    var y1r = .0
    var y1i = .0
    var y2r = .0
    var y2i = .0
    var y3r = .0
    var y3i = .0
    var y4r = .0
    var y4i = .0
    var y5r = .0
    var y5i = .0
    var y6r = .0
    var y6i = .0
    var y7r = .0
    var y7i = .0
    var y8r = .0
    var y8i = .0
    var y9r = .0
    var y9i = .0
    var y10r = .0
    var y10i = .0
    var y11r = .0
    var y11i = .0
    var y12r = .0
    var y12i = .0
    var y13r = .0
    var y13i = .0
    var y14r = .0
    var y14i = .0
    var y15r = .0
    var y15i = .0
    wn4r = w(startw + 1)
    wk1r = w(startw + 4)
    wk1i = w(startw + 5)
    wk3r = w(startw + 6)
    wk3i = -w(startw + 7)
    wk2r = w(startw + 8)
    wk2i = w(startw + 9)
    x1r = a(offa) - a(offa + 17)
    x1i = a(offa + 1) + a(offa + 16)
    x0r = a(offa + 8) - a(offa + 25)
    x0i = a(offa + 9) + a(offa + 24)
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    y0r = x1r + x2r
    y0i = x1i + x2i
    y4r = x1r - x2r
    y4i = x1i - x2i
    x1r = a(offa) + a(offa + 17)
    x1i = a(offa + 1) - a(offa + 16)
    x0r = a(offa + 8) + a(offa + 25)
    x0i = a(offa + 9) - a(offa + 24)
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    y8r = x1r - x2i
    y8i = x1i + x2r
    y12r = x1r + x2i
    y12i = x1i - x2r
    x0r = a(offa + 2) - a(offa + 19)
    x0i = a(offa + 3) + a(offa + 18)
    x1r = wk1r * x0r - wk1i * x0i
    x1i = wk1r * x0i + wk1i * x0r
    x0r = a(offa + 10) - a(offa + 27)
    x0i = a(offa + 11) + a(offa + 26)
    x2r = wk3i * x0r - wk3r * x0i
    x2i = wk3i * x0i + wk3r * x0r
    y1r = x1r + x2r
    y1i = x1i + x2i
    y5r = x1r - x2r
    y5i = x1i - x2i
    x0r = a(offa + 2) + a(offa + 19)
    x0i = a(offa + 3) - a(offa + 18)
    x1r = wk3r * x0r - wk3i * x0i
    x1i = wk3r * x0i + wk3i * x0r
    x0r = a(offa + 10) + a(offa + 27)
    x0i = a(offa + 11) - a(offa + 26)
    x2r = wk1r * x0r + wk1i * x0i
    x2i = wk1r * x0i - wk1i * x0r
    y9r = x1r - x2r
    y9i = x1i - x2i
    y13r = x1r + x2r
    y13i = x1i + x2i
    x0r = a(offa + 4) - a(offa + 21)
    x0i = a(offa + 5) + a(offa + 20)
    x1r = wk2r * x0r - wk2i * x0i
    x1i = wk2r * x0i + wk2i * x0r
    x0r = a(offa + 12) - a(offa + 29)
    x0i = a(offa + 13) + a(offa + 28)
    x2r = wk2i * x0r - wk2r * x0i
    x2i = wk2i * x0i + wk2r * x0r
    y2r = x1r + x2r
    y2i = x1i + x2i
    y6r = x1r - x2r
    y6i = x1i - x2i
    x0r = a(offa + 4) + a(offa + 21)
    x0i = a(offa + 5) - a(offa + 20)
    x1r = wk2i * x0r - wk2r * x0i
    x1i = wk2i * x0i + wk2r * x0r
    x0r = a(offa + 12) + a(offa + 29)
    x0i = a(offa + 13) - a(offa + 28)
    x2r = wk2r * x0r - wk2i * x0i
    x2i = wk2r * x0i + wk2i * x0r
    y10r = x1r - x2r
    y10i = x1i - x2i
    y14r = x1r + x2r
    y14i = x1i + x2i
    x0r = a(offa + 6) - a(offa + 23)
    x0i = a(offa + 7) + a(offa + 22)
    x1r = wk3r * x0r - wk3i * x0i
    x1i = wk3r * x0i + wk3i * x0r
    x0r = a(offa + 14) - a(offa + 31)
    x0i = a(offa + 15) + a(offa + 30)
    x2r = wk1i * x0r - wk1r * x0i
    x2i = wk1i * x0i + wk1r * x0r
    y3r = x1r + x2r
    y3i = x1i + x2i
    y7r = x1r - x2r
    y7i = x1i - x2i
    x0r = a(offa + 6) + a(offa + 23)
    x0i = a(offa + 7) - a(offa + 22)
    x1r = wk1i * x0r + wk1r * x0i
    x1i = wk1i * x0i - wk1r * x0r
    x0r = a(offa + 14) + a(offa + 31)
    x0i = a(offa + 15) - a(offa + 30)
    x2r = wk3i * x0r - wk3r * x0i
    x2i = wk3i * x0i + wk3r * x0r
    y11r = x1r + x2r
    y11i = x1i + x2i
    y15r = x1r - x2r
    y15i = x1i - x2i
    x1r = y0r + y2r
    x1i = y0i + y2i
    x2r = y1r + y3r
    x2i = y1i + y3i
    a(offa) = x1r + x2r
    a(offa + 1) = x1i + x2i
    a(offa + 2) = x1r - x2r
    a(offa + 3) = x1i - x2i
    x1r = y0r - y2r
    x1i = y0i - y2i
    x2r = y1r - y3r
    x2i = y1i - y3i
    a(offa + 4) = x1r - x2i
    a(offa + 5) = x1i + x2r
    a(offa + 6) = x1r + x2i
    a(offa + 7) = x1i - x2r
    x1r = y4r - y6i
    x1i = y4i + y6r
    x0r = y5r - y7i
    x0i = y5i + y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 8) = x1r + x2r
    a(offa + 9) = x1i + x2i
    a(offa + 10) = x1r - x2r
    a(offa + 11) = x1i - x2i
    x1r = y4r + y6i
    x1i = y4i - y6r
    x0r = y5r + y7i
    x0i = y5i - y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 12) = x1r - x2i
    a(offa + 13) = x1i + x2r
    a(offa + 14) = x1r + x2i
    a(offa + 15) = x1i - x2r
    x1r = y8r + y10r
    x1i = y8i + y10i
    x2r = y9r - y11r
    x2i = y9i - y11i
    a(offa + 16) = x1r + x2r
    a(offa + 17) = x1i + x2i
    a(offa + 18) = x1r - x2r
    a(offa + 19) = x1i - x2i
    x1r = y8r - y10r
    x1i = y8i - y10i
    x2r = y9r + y11r
    x2i = y9i + y11i
    a(offa + 20) = x1r - x2i
    a(offa + 21) = x1i + x2r
    a(offa + 22) = x1r + x2i
    a(offa + 23) = x1i - x2r
    x1r = y12r - y14i
    x1i = y12i + y14r
    x0r = y13r + y15i
    x0i = y13i - y15r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 24) = x1r + x2r
    a(offa + 25) = x1i + x2i
    a(offa + 26) = x1r - x2r
    a(offa + 27) = x1i - x2i
    x1r = y12r + y14i
    x1i = y12i - y14r
    x0r = y13r - y15i
    x0i = y13i + y15r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 28) = x1r - x2i
    a(offa + 29) = x1i + x2r
    a(offa + 30) = x1r + x2i
    a(offa + 31) = x1i - x2r
  }

  def cftf081(a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var wn4r = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    var y0r = .0
    var y0i = .0
    var y1r = .0
    var y1i = .0
    var y2r = .0
    var y2i = .0
    var y3r = .0
    var y3i = .0
    var y4r = .0
    var y4i = .0
    var y5r = .0
    var y5i = .0
    var y6r = .0
    var y6i = .0
    var y7r = .0
    var y7i = .0
    wn4r = w(startw + 1)
    x0r = a(offa) + a(offa + 8)
    x0i = a(offa + 1) + a(offa + 9)
    x1r = a(offa) - a(offa + 8)
    x1i = a(offa + 1) - a(offa + 9)
    x2r = a(offa + 4) + a(offa + 12)
    x2i = a(offa + 5) + a(offa + 13)
    x3r = a(offa + 4) - a(offa + 12)
    x3i = a(offa + 5) - a(offa + 13)
    y0r = x0r + x2r
    y0i = x0i + x2i
    y2r = x0r - x2r
    y2i = x0i - x2i
    y1r = x1r - x3i
    y1i = x1i + x3r
    y3r = x1r + x3i
    y3i = x1i - x3r
    x0r = a(offa + 2) + a(offa + 10)
    x0i = a(offa + 3) + a(offa + 11)
    x1r = a(offa + 2) - a(offa + 10)
    x1i = a(offa + 3) - a(offa + 11)
    x2r = a(offa + 6) + a(offa + 14)
    x2i = a(offa + 7) + a(offa + 15)
    x3r = a(offa + 6) - a(offa + 14)
    x3i = a(offa + 7) - a(offa + 15)
    y4r = x0r + x2r
    y4i = x0i + x2i
    y6r = x0r - x2r
    y6i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    x2r = x1r + x3i
    x2i = x1i - x3r
    y5r = wn4r * (x0r - x0i)
    y5i = wn4r * (x0r + x0i)
    y7r = wn4r * (x2r - x2i)
    y7i = wn4r * (x2r + x2i)
    a(offa + 8) = y1r + y5r
    a(offa + 9) = y1i + y5i
    a(offa + 10) = y1r - y5r
    a(offa + 11) = y1i - y5i
    a(offa + 12) = y3r - y7i
    a(offa + 13) = y3i + y7r
    a(offa + 14) = y3r + y7i
    a(offa + 15) = y3i - y7r
    a(offa) = y0r + y4r
    a(offa + 1) = y0i + y4i
    a(offa + 2) = y0r - y4r
    a(offa + 3) = y0i - y4i
    a(offa + 4) = y2r - y6i
    a(offa + 5) = y2i + y6r
    a(offa + 6) = y2r + y6i
    a(offa + 7) = y2i - y6r
  }

  def cftf082(a: Array[Double], offa: Int, w: Array[Double], startw: Int): Unit = {
    var wn4r = .0
    var wk1r = .0
    var wk1i = .0
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var y0r = .0
    var y0i = .0
    var y1r = .0
    var y1i = .0
    var y2r = .0
    var y2i = .0
    var y3r = .0
    var y3i = .0
    var y4r = .0
    var y4i = .0
    var y5r = .0
    var y5i = .0
    var y6r = .0
    var y6i = .0
    var y7r = .0
    var y7i = .0
    wn4r = w(startw + 1)
    wk1r = w(startw + 2)
    wk1i = w(startw + 3)
    y0r = a(offa) - a(offa + 9)
    y0i = a(offa + 1) + a(offa + 8)
    y1r = a(offa) + a(offa + 9)
    y1i = a(offa + 1) - a(offa + 8)
    x0r = a(offa + 4) - a(offa + 13)
    x0i = a(offa + 5) + a(offa + 12)
    y2r = wn4r * (x0r - x0i)
    y2i = wn4r * (x0i + x0r)
    x0r = a(offa + 4) + a(offa + 13)
    x0i = a(offa + 5) - a(offa + 12)
    y3r = wn4r * (x0r - x0i)
    y3i = wn4r * (x0i + x0r)
    x0r = a(offa + 2) - a(offa + 11)
    x0i = a(offa + 3) + a(offa + 10)
    y4r = wk1r * x0r - wk1i * x0i
    y4i = wk1r * x0i + wk1i * x0r
    x0r = a(offa + 2) + a(offa + 11)
    x0i = a(offa + 3) - a(offa + 10)
    y5r = wk1i * x0r - wk1r * x0i
    y5i = wk1i * x0i + wk1r * x0r
    x0r = a(offa + 6) - a(offa + 15)
    x0i = a(offa + 7) + a(offa + 14)
    y6r = wk1i * x0r - wk1r * x0i
    y6i = wk1i * x0i + wk1r * x0r
    x0r = a(offa + 6) + a(offa + 15)
    x0i = a(offa + 7) - a(offa + 14)
    y7r = wk1r * x0r - wk1i * x0i
    y7i = wk1r * x0i + wk1i * x0r
    x0r = y0r + y2r
    x0i = y0i + y2i
    x1r = y4r + y6r
    x1i = y4i + y6i
    a(offa) = x0r + x1r
    a(offa + 1) = x0i + x1i
    a(offa + 2) = x0r - x1r
    a(offa + 3) = x0i - x1i
    x0r = y0r - y2r
    x0i = y0i - y2i
    x1r = y4r - y6r
    x1i = y4i - y6i
    a(offa + 4) = x0r - x1i
    a(offa + 5) = x0i + x1r
    a(offa + 6) = x0r + x1i
    a(offa + 7) = x0i - x1r
    x0r = y1r - y3i
    x0i = y1i + y3r
    x1r = y5r - y7r
    x1i = y5i - y7i
    a(offa + 8) = x0r + x1r
    a(offa + 9) = x0i + x1i
    a(offa + 10) = x0r - x1r
    a(offa + 11) = x0i - x1i
    x0r = y1r + y3i
    x0i = y1i - y3r
    x1r = y5r + y7r
    x1i = y5i + y7i
    a(offa + 12) = x0r - x1i
    a(offa + 13) = x0i + x1r
    a(offa + 14) = x0r + x1i
    a(offa + 15) = x0i - x1r
  }

  def cftf040(a: Array[Double], offa: Int): Unit = {
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    x0r = a(offa) + a(offa + 4)
    x0i = a(offa + 1) + a(offa + 5)
    x1r = a(offa) - a(offa + 4)
    x1i = a(offa + 1) - a(offa + 5)
    x2r = a(offa + 2) + a(offa + 6)
    x2i = a(offa + 3) + a(offa + 7)
    x3r = a(offa + 2) - a(offa + 6)
    x3i = a(offa + 3) - a(offa + 7)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(offa + 2) = x1r - x3i
    a(offa + 3) = x1i + x3r
    a(offa + 4) = x0r - x2r
    a(offa + 5) = x0i - x2i
    a(offa + 6) = x1r + x3i
    a(offa + 7) = x1i - x3r
  }

  def cftb040(a: Array[Double], offa: Int): Unit = {
    var x0r = .0
    var x0i = .0
    var x1r = .0
    var x1i = .0
    var x2r = .0
    var x2i = .0
    var x3r = .0
    var x3i = .0
    x0r = a(offa) + a(offa + 4)
    x0i = a(offa + 1) + a(offa + 5)
    x1r = a(offa) - a(offa + 4)
    x1i = a(offa + 1) - a(offa + 5)
    x2r = a(offa + 2) + a(offa + 6)
    x2i = a(offa + 3) + a(offa + 7)
    x3r = a(offa + 2) - a(offa + 6)
    x3i = a(offa + 3) - a(offa + 7)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(offa + 2) = x1r + x3i
    a(offa + 3) = x1i - x3r
    a(offa + 4) = x0r - x2r
    a(offa + 5) = x0i - x2i
    a(offa + 6) = x1r - x3i
    a(offa + 7) = x1i + x3r
  }

  def cftx020(a: Array[Double], offa: Int): Unit = {
    var x0r = .0
    var x0i = .0
    x0r = a(offa) - a(offa + 2)
    x0i = -a(offa + 1) + a(offa + 3)
    a(offa) += a(offa + 2)
    a(offa + 1) += a(offa + 3)
    a(offa + 2) = x0r
    a(offa + 3) = x0i
  }

  def cftxb020(a: Array[Double], offa: Int): Unit = {
    var x0r = .0
    var x0i = .0
    x0r = a(offa) - a(offa + 2)
    x0i = a(offa + 1) - a(offa + 3)
    a(offa) += a(offa + 2)
    a(offa + 1) += a(offa + 3)
    a(offa + 2) = x0r
    a(offa + 3) = x0i
  }

  def cftxc020(a: Array[Double], offa: Int): Unit = {
    var x0r = .0
    var x0i = .0
    x0r = a(offa) - a(offa + 2)
    x0i = a(offa + 1) + a(offa + 3)
    a(offa) += a(offa + 2)
    a(offa + 1) -= a(offa + 3)
    a(offa + 2) = x0r
    a(offa + 3) = x0i
  }

  def rftfsub(n: Int, a: Array[Double], offa: Int, nc: Int, c: Array[Double], startc: Int): Unit = {
    var k = 0
    var kk = 0
    var ks = 0
    var m = 0
    var wkr = .0
    var wki = .0
    var xr = .0
    var xi = .0
    var yr = .0
    var yi = .0
    var idx1 = 0
    var idx2 = 0
    m = n >> 1
    ks = 2 * nc / m
    kk = 0
    var j = 2
    while ( {
      j < m
    }) {
      k = n - j
      kk += ks
      wkr = 0.5 - c(startc + nc - kk)
      wki = c(startc + kk)
      idx1 = offa + j
      idx2 = offa + k
      xr = a(idx1) - a(idx2)
      xi = a(idx1 + 1) + a(idx2 + 1)
      yr = wkr * xr - wki * xi
      yi = wkr * xi + wki * xr
      a(idx1) -= yr
      a(idx1 + 1) = yi - a(idx1 + 1)
      a(idx2) += yr
      a(idx2 + 1) = yi - a(idx2 + 1)

      j += 2
    }
    a(offa + m + 1) = -a(offa + m + 1)
  }

  def rftbsub(n: Int, a: Array[Double], offa: Int, nc: Int, c: Array[Double], startc: Int): Unit = {
    var k = 0
    var kk = 0
    var ks = 0
    var m = 0
    var wkr = .0
    var wki = .0
    var xr = .0
    var xi = .0
    var yr = .0
    var yi = .0
    var idx1 = 0
    var idx2 = 0
    m = n >> 1
    ks = 2 * nc / m
    kk = 0
    var j = 2
    while ( {
      j < m
    }) {
      k = n - j
      kk += ks
      wkr = 0.5 - c(startc + nc - kk)
      wki = c(startc + kk)
      idx1 = offa + j
      idx2 = offa + k
      xr = a(idx1) - a(idx2)
      xi = a(idx1 + 1) + a(idx2 + 1)
      yr = wkr * xr - wki * xi
      yi = wkr * xi + wki * xr
      a(idx1) -= yr
      a(idx1 + 1) -= yi
      a(idx2) += yr
      a(idx2 + 1) -= yi

      j += 2
    }
  }

  def dctsub(n: Int, a: Array[Double], offa: Int, nc: Int, c: Array[Double], startc: Int): Unit = {
    var k = 0
    var kk = 0
    var ks = 0
    var m = 0
    var wkr = .0
    var wki = .0
    var xr = .0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    m = n >> 1
    ks = nc / n
    kk = 0
    for (j <- 1 until m) {
      k = n - j
      kk += ks
      idx0 = startc + kk
      idx1 = offa + j
      idx2 = offa + k
      wkr = c(idx0) - c(startc + nc - kk)
      wki = c(idx0) + c(startc + nc - kk)
      xr = wki * a(idx1) - wkr * a(idx2)
      a(idx1) = wkr * a(idx1) + wki * a(idx2)
      a(idx2) = xr
    }
    a(offa + m) *= c(startc)
  }

  def cftfsub(n: Int, a: Array[Float], offa: Int, ip: Array[Int], nw: Int, w: Array[Float]): Unit = {
    if (n > 8) if (n > 32) {
      cftf1st(n, a, offa, w, nw - (n >> 2))
      if ((ConcurrencyUtils.numThreads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) cftrec4_th(n, a, offa, nw, w)
      else if (n > 512) cftrec4(n, a, offa, nw, w)
      else if (n > 128) cftleaf(n, 1, a, offa, nw, w)
      else cftfx41(n, a, offa, nw, w)
      bitrv2(n, ip, a, offa)
    }
    else if (n == 32) {
      cftf161(a, offa, w, nw - 8)
      bitrv216(a, offa)
    }
    else {
      cftf081(a, offa, w, 0)
      bitrv208(a, offa)
    }
    else if (n == 8) cftf040(a, offa)
    else if (n == 4) cftxb020(a, offa)
  }

  def cftbsub(n: Int, a: Array[Float], offa: Int, ip: Array[Int], nw: Int, w: Array[Float]): Unit = {
    if (n > 8) if (n > 32) {
      cftb1st(n, a, offa, w, nw - (n >> 2))
      if ((ConcurrencyUtils.numThreads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) cftrec4_th(n, a, offa, nw, w)
      else if (n > 512) cftrec4(n, a, offa, nw, w)
      else if (n > 128) cftleaf(n, 1, a, offa, nw, w)
      else cftfx41(n, a, offa, nw, w)
      bitrv2conj(n, ip, a, offa)
    }
    else if (n == 32) {
      cftf161(a, offa, w, nw - 8)
      bitrv216neg(a, offa)
    }
    else {
      cftf081(a, offa, w, 0)
      bitrv208neg(a, offa)
    }
    else if (n == 8) cftb040(a, offa)
    else if (n == 4) cftxb020(a, offa)
  }

  def bitrv2(n: Int, ip: Array[Int], a: Array[Float], offa: Int): Unit = {
    var j1    = 0
    var k1    = 0
    var l     = 0
    var m     = 0
    var nh    = 0
    var nm    = 0
    var xr    = 0.0f
    var xi    = 0.0f
    var yr    = 0.0f
    var yi    = 0.0f
    var idx0  = 0
    var idx1  = 0
    var idx2  = 0
    m = 1
    l = n >> 2
    while ( {
      l > 8
    }) {
      m <<= 1

      l >>= 2
    }
    nh = n >> 1
    nm = 4 * m
    if (l == 8) for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + 2 * ip(m + k)
        k1 = idx0 + 2 * ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + 2 * ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 += 2 * nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 -= nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= 2
      k1 -= nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nh + 2
      k1 += nh + 2
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= nh - nm
      k1 += 2 * nm - 2
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
    }
    else for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + ip(m + k)
        k1 = idx0 + ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = a(idx1 + 1)
        yr = a(idx2)
        yi = a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 += nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = a(idx1 + 1)
      yr = a(idx2)
      yi = a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
    }
  }

  def bitrv2conj(n: Int, ip: Array[Int], a: Array[Float], offa: Int): Unit = {
    var j1    = 0
    var k1    = 0
    var l     = 0
    var m     = 0
    var nh    = 0
    var nm    = 0
    var xr    = 0.0f
    var xi    = 0.0f
    var yr    = 0.0f
    var yi    = 0.0f
    var idx0  = 0
    var idx1  = 0
    var idx2  = 0
    m = 1
    l = n >> 2
    while ( {
      l > 8
    }) {
      m <<= 1

      l >>= 2
    }
    nh = n >> 1
    nm = 4 * m
    if (l == 8) for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + 2 * ip(m + k)
        k1 = idx0 + 2 * ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= 2 * nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + 2 * ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
      j1 += nm
      k1 += 2 * nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nm
      k1 -= nm
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= 2
      k1 -= nh
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 += nh + 2
      k1 += nh + 2
      idx1 = offa + j1
      idx2 = offa + k1
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      j1 -= nh - nm
      k1 += 2 * nm - 2
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
    }
    else for (k <- 0 until m) {
      idx0 = 4 * k
      for (j <- 0 until k) {
        j1 = 4 * j + ip(m + k)
        k1 = idx0 + ip(m + j)
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nh
        k1 += 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += 2
        k1 += nh
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 += nm
        k1 += nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nh
        k1 -= 2
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
        j1 -= nm
        k1 -= nm
        idx1 = offa + j1
        idx2 = offa + k1
        xr = a(idx1)
        xi = -a(idx1 + 1)
        yr = a(idx2)
        yi = -a(idx2 + 1)
        a(idx1) = yr
        a(idx1 + 1) = yi
        a(idx2) = xr
        a(idx2 + 1) = xi
      }
      k1 = idx0 + ip(m + k)
      j1 = k1 + 2
      k1 += nh
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
      j1 += nm
      k1 += nm
      idx1 = offa + j1
      idx2 = offa + k1
      a(idx1 - 1) = -a(idx1 - 1)
      xr = a(idx1)
      xi = -a(idx1 + 1)
      yr = a(idx2)
      yi = -a(idx2 + 1)
      a(idx1) = yr
      a(idx1 + 1) = yi
      a(idx2) = xr
      a(idx2 + 1) = xi
      a(idx2 + 3) = -a(idx2 + 3)
    }
  }

  def bitrv216(a: Array[Float], offa: Int): Unit = {
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var x4r   = 0.0f
    var x4i   = 0.0f
    var x5r   = 0.0f
    var x5i   = 0.0f
    var x7r   = 0.0f
    var x7i   = 0.0f
    var x8r   = 0.0f
    var x8i   = 0.0f
    var x10r  = 0.0f
    var x10i  = 0.0f
    var x11r  = 0.0f
    var x11i  = 0.0f
    var x12r  = 0.0f
    var x12i  = 0.0f
    var x13r  = 0.0f
    var x13i  = 0.0f
    var x14r  = 0.0f
    var x14i  = 0.0f
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x2r = a(offa + 4)
    x2i = a(offa + 5)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x5r = a(offa + 10)
    x5i = a(offa + 11)
    x7r = a(offa + 14)
    x7i = a(offa + 15)
    x8r = a(offa + 16)
    x8i = a(offa + 17)
    x10r = a(offa + 20)
    x10i = a(offa + 21)
    x11r = a(offa + 22)
    x11i = a(offa + 23)
    x12r = a(offa + 24)
    x12i = a(offa + 25)
    x13r = a(offa + 26)
    x13i = a(offa + 27)
    x14r = a(offa + 28)
    x14i = a(offa + 29)
    a(offa + 2) = x8r
    a(offa + 3) = x8i
    a(offa + 4) = x4r
    a(offa + 5) = x4i
    a(offa + 6) = x12r
    a(offa + 7) = x12i
    a(offa + 8) = x2r
    a(offa + 9) = x2i
    a(offa + 10) = x10r
    a(offa + 11) = x10i
    a(offa + 14) = x14r
    a(offa + 15) = x14i
    a(offa + 16) = x1r
    a(offa + 17) = x1i
    a(offa + 20) = x5r
    a(offa + 21) = x5i
    a(offa + 22) = x13r
    a(offa + 23) = x13i
    a(offa + 24) = x3r
    a(offa + 25) = x3i
    a(offa + 26) = x11r
    a(offa + 27) = x11i
    a(offa + 28) = x7r
    a(offa + 29) = x7i
  }

  def bitrv216neg(a: Array[Float], offa: Int): Unit = {
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var x4r   = 0.0f
    var x4i   = 0.0f
    var x5r   = 0.0f
    var x5i   = 0.0f
    var x6r   = 0.0f
    var x6i   = 0.0f
    var x7r   = 0.0f
    var x7i   = 0.0f
    var x8r   = 0.0f
    var x8i   = 0.0f
    var x9r   = 0.0f
    var x9i   = 0.0f
    var x10r  = 0.0f
    var x10i  = 0.0f
    var x11r  = 0.0f
    var x11i  = 0.0f
    var x12r  = 0.0f
    var x12i  = 0.0f
    var x13r  = 0.0f
    var x13i  = 0.0f
    var x14r  = 0.0f
    var x14i  = 0.0f
    var x15r  = 0.0f
    var x15i  = 0.0f
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x2r = a(offa + 4)
    x2i = a(offa + 5)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x5r = a(offa + 10)
    x5i = a(offa + 11)
    x6r = a(offa + 12)
    x6i = a(offa + 13)
    x7r = a(offa + 14)
    x7i = a(offa + 15)
    x8r = a(offa + 16)
    x8i = a(offa + 17)
    x9r = a(offa + 18)
    x9i = a(offa + 19)
    x10r = a(offa + 20)
    x10i = a(offa + 21)
    x11r = a(offa + 22)
    x11i = a(offa + 23)
    x12r = a(offa + 24)
    x12i = a(offa + 25)
    x13r = a(offa + 26)
    x13i = a(offa + 27)
    x14r = a(offa + 28)
    x14i = a(offa + 29)
    x15r = a(offa + 30)
    x15i = a(offa + 31)
    a(offa + 2) = x15r
    a(offa + 3) = x15i
    a(offa + 4) = x7r
    a(offa + 5) = x7i
    a(offa + 6) = x11r
    a(offa + 7) = x11i
    a(offa + 8) = x3r
    a(offa + 9) = x3i
    a(offa + 10) = x13r
    a(offa + 11) = x13i
    a(offa + 12) = x5r
    a(offa + 13) = x5i
    a(offa + 14) = x9r
    a(offa + 15) = x9i
    a(offa + 16) = x1r
    a(offa + 17) = x1i
    a(offa + 18) = x14r
    a(offa + 19) = x14i
    a(offa + 20) = x6r
    a(offa + 21) = x6i
    a(offa + 22) = x10r
    a(offa + 23) = x10i
    a(offa + 24) = x2r
    a(offa + 25) = x2i
    a(offa + 26) = x12r
    a(offa + 27) = x12i
    a(offa + 28) = x4r
    a(offa + 29) = x4i
    a(offa + 30) = x8r
    a(offa + 31) = x8i
  }

  def bitrv208(a: Array[Float], offa: Int): Unit = {
    var x1r = 0.0f
    var x1i = 0.0f
    var x3r = 0.0f
    var x3i = 0.0f
    var x4r = 0.0f
    var x4i = 0.0f
    var x6r = 0.0f
    var x6i = 0.0f
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x6r = a(offa + 12)
    x6i = a(offa + 13)
    a(offa + 2) = x4r
    a(offa + 3) = x4i
    a(offa + 6) = x6r
    a(offa + 7) = x6i
    a(offa + 8) = x1r
    a(offa + 9) = x1i
    a(offa + 12) = x3r
    a(offa + 13) = x3i
  }

  def bitrv208neg(a: Array[Float], offa: Int): Unit = {
    var x1r = 0.0f
    var x1i = 0.0f
    var x2r = 0.0f
    var x2i = 0.0f
    var x3r = 0.0f
    var x3i = 0.0f
    var x4r = 0.0f
    var x4i = 0.0f
    var x5r = 0.0f
    var x5i = 0.0f
    var x6r = 0.0f
    var x6i = 0.0f
    var x7r = 0.0f
    var x7i = 0.0f
    x1r = a(offa + 2)
    x1i = a(offa + 3)
    x2r = a(offa + 4)
    x2i = a(offa + 5)
    x3r = a(offa + 6)
    x3i = a(offa + 7)
    x4r = a(offa + 8)
    x4i = a(offa + 9)
    x5r = a(offa + 10)
    x5i = a(offa + 11)
    x6r = a(offa + 12)
    x6i = a(offa + 13)
    x7r = a(offa + 14)
    x7i = a(offa + 15)
    a(offa + 2) = x7r
    a(offa + 3) = x7i
    a(offa + 4) = x3r
    a(offa + 5) = x3i
    a(offa + 6) = x5r
    a(offa + 7) = x5i
    a(offa + 8) = x1r
    a(offa + 9) = x1i
    a(offa + 10) = x6r
    a(offa + 11) = x6i
    a(offa + 12) = x2r
    a(offa + 13) = x2i
    a(offa + 14) = x4r
    a(offa + 15) = x4i
  }

  def cftf1st(n: Int, a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var j0    = 0
    var j1    = 0
    var j2    = 0
    var j3    = 0
    var k     = 0
    var m     = 0
    var mh    = 0
    var wn4r  = 0.0f
    var csc1  = 0.0f
    var csc3  = 0.0f
    var wk1r  = 0.0f
    var wk1i  = 0.0f
    var wk3r  = 0.0f
    var wk3i  = 0.0f
    var wd1r  = 0.0f
    var wd1i  = 0.0f
    var wd3r  = 0.0f
    var wd3i  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var y0r   = 0.0f
    var y0i   = 0.0f
    var y1r   = 0.0f
    var y1i   = 0.0f
    var y2r   = 0.0f
    var y2i   = 0.0f
    var y3r   = 0.0f
    var y3i   = 0.0f
    var idx0  = 0
    var idx1  = 0
    var idx2  = 0
    var idx3  = 0
    var idx4  = 0
    var idx5  = 0
    mh = n >> 3
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) + a(idx2)
    x0i = a(offa + 1) + a(idx2 + 1)
    x1r = a(offa) - a(idx2)
    x1i = a(offa + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    a(idx2) = x1r - x3i
    a(idx2 + 1) = x1i + x3r
    a(idx3) = x1r + x3i
    a(idx3 + 1) = x1i - x3r
    wn4r = w(startw + 1)
    csc1 = w(startw + 2)
    csc3 = w(startw + 3)
    wd1r = 1
    wd1i = 0
    wd3r = 1
    wd3i = 0
    k = 0
    var j = 2
    while ( {
      j < mh - 2
    }) {
      k += 4
      idx4 = startw + k
      wk1r = csc1 * (wd1r + w(idx4))
      wk1i = csc1 * (wd1i + w(idx4 + 1))
      wk3r = csc3 * (wd3r + w(idx4 + 2))
      wk3i = csc3 * (wd3i + w(idx4 + 3))
      wd1r = w(idx4)
      wd1i = w(idx4 + 1)
      wd3r = w(idx4 + 2)
      wd3i = w(idx4 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx5 = offa + j
      x0r = a(idx5) + a(idx2)
      x0i = a(idx5 + 1) + a(idx2 + 1)
      x1r = a(idx5) - a(idx2)
      x1i = a(idx5 + 1) - a(idx2 + 1)
      y0r = a(idx5 + 2) + a(idx2 + 2)
      y0i = a(idx5 + 3) + a(idx2 + 3)
      y1r = a(idx5 + 2) - a(idx2 + 2)
      y1i = a(idx5 + 3) - a(idx2 + 3)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 + 2) + a(idx3 + 2)
      y2i = a(idx1 + 3) + a(idx3 + 3)
      y3r = a(idx1 + 2) - a(idx3 + 2)
      y3i = a(idx1 + 3) - a(idx3 + 3)
      a(idx5) = x0r + x2r
      a(idx5 + 1) = x0i + x2i
      a(idx5 + 2) = y0r + y2r
      a(idx5 + 3) = y0i + y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      a(idx1 + 2) = y0r - y2r
      a(idx1 + 3) = y0i - y2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1r * x0r - wk1i * x0i
      a(idx2 + 1) = wk1r * x0i + wk1i * x0r
      x0r = y1r - y3i
      x0i = y1i + y3r
      a(idx2 + 2) = wd1r * x0r - wd1i * x0i
      a(idx2 + 3) = wd1r * x0i + wd1i * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3r * x0r + wk3i * x0i
      a(idx3 + 1) = wk3r * x0i - wk3i * x0r
      x0r = y1r + y3i
      x0i = y1i - y3r
      a(idx3 + 2) = wd3r * x0r + wd3i * x0i
      a(idx3 + 3) = wd3r * x0i - wd3i * x0r
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) + a(idx2)
      x0i = a(idx0 + 1) + a(idx2 + 1)
      x1r = a(idx0) - a(idx2)
      x1i = a(idx0 + 1) - a(idx2 + 1)
      y0r = a(idx0 - 2) + a(idx2 - 2)
      y0i = a(idx0 - 1) + a(idx2 - 1)
      y1r = a(idx0 - 2) - a(idx2 - 2)
      y1i = a(idx0 - 1) - a(idx2 - 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 - 2) + a(idx3 - 2)
      y2i = a(idx1 - 1) + a(idx3 - 1)
      y3r = a(idx1 - 2) - a(idx3 - 2)
      y3i = a(idx1 - 1) - a(idx3 - 1)
      a(idx0) = x0r + x2r
      a(idx0 + 1) = x0i + x2i
      a(idx0 - 2) = y0r + y2r
      a(idx0 - 1) = y0i + y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      a(idx1 - 2) = y0r - y2r
      a(idx1 - 1) = y0i - y2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1i * x0r - wk1r * x0i
      a(idx2 + 1) = wk1i * x0i + wk1r * x0r
      x0r = y1r - y3i
      x0i = y1i + y3r
      a(idx2 - 2) = wd1i * x0r - wd1r * x0i
      a(idx2 - 1) = wd1i * x0i + wd1r * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3i * x0r + wk3r * x0i
      a(idx3 + 1) = wk3i * x0i - wk3r * x0r
      x0r = y1r + y3i
      x0i = y1i - y3r
      a(offa + j3 - 2) = wd3i * x0r + wd3r * x0i
      a(offa + j3 - 1) = wd3i * x0i - wd3r * x0r

      j += 4
    }
    wk1r = csc1 * (wd1r + wn4r)
    wk1i = csc1 * (wd1i + wn4r)
    wk3r = csc3 * (wd3r - wn4r)
    wk3i = csc3 * (wd3i - wn4r)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0 - 2) + a(idx2 - 2)
    x0i = a(idx0 - 1) + a(idx2 - 1)
    x1r = a(idx0 - 2) - a(idx2 - 2)
    x1i = a(idx0 - 1) - a(idx2 - 1)
    x2r = a(idx1 - 2) + a(idx3 - 2)
    x2i = a(idx1 - 1) + a(idx3 - 1)
    x3r = a(idx1 - 2) - a(idx3 - 2)
    x3i = a(idx1 - 1) - a(idx3 - 1)
    a(idx0 - 2) = x0r + x2r
    a(idx0 - 1) = x0i + x2i
    a(idx1 - 2) = x0r - x2r
    a(idx1 - 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2 - 2) = wk1r * x0r - wk1i * x0i
    a(idx2 - 1) = wk1r * x0i + wk1i * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3 - 2) = wk3r * x0r + wk3i * x0i
    a(idx3 - 1) = wk3r * x0i - wk3i * x0r
    x0r = a(idx0) + a(idx2)
    x0i = a(idx0 + 1) + a(idx2 + 1)
    x1r = a(idx0) - a(idx2)
    x1i = a(idx0 + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(idx0) = x0r + x2r
    a(idx0 + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2) = wn4r * (x0r - x0i)
    a(idx2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3) = -wn4r * (x0r + x0i)
    a(idx3 + 1) = -wn4r * (x0i - x0r)
    x0r = a(idx0 + 2) + a(idx2 + 2)
    x0i = a(idx0 + 3) + a(idx2 + 3)
    x1r = a(idx0 + 2) - a(idx2 + 2)
    x1i = a(idx0 + 3) - a(idx2 + 3)
    x2r = a(idx1 + 2) + a(idx3 + 2)
    x2i = a(idx1 + 3) + a(idx3 + 3)
    x3r = a(idx1 + 2) - a(idx3 + 2)
    x3i = a(idx1 + 3) - a(idx3 + 3)
    a(idx0 + 2) = x0r + x2r
    a(idx0 + 3) = x0i + x2i
    a(idx1 + 2) = x0r - x2r
    a(idx1 + 3) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2 + 2) = wk1i * x0r - wk1r * x0i
    a(idx2 + 3) = wk1i * x0i + wk1r * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3 + 2) = wk3i * x0r + wk3r * x0i
    a(idx3 + 3) = wk3i * x0i - wk3r * x0r
  }

  def cftb1st(n: Int, a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var j0    = 0
    var j1    = 0
    var j2    = 0
    var j3    = 0
    var k     = 0
    var m     = 0
    var mh    = 0
    var wn4r  = 0.0f
    var csc1  = 0.0f
    var csc3  = 0.0f
    var wk1r  = 0.0f
    var wk1i  = 0.0f
    var wk3r  = 0.0f
    var wk3i  = 0.0f
    var wd1r  = 0.0f
    var wd1i  = 0.0f
    var wd3r  = 0.0f
    var wd3i  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var y0r   = 0.0f
    var y0i   = 0.0f
    var y1r   = 0.0f
    var y1i   = 0.0f
    var y2r   = 0.0f
    var y2i   = 0.0f
    var y3r   = 0.0f
    var y3i   = 0.0f
    var idx0  = 0
    var idx1  = 0
    var idx2  = 0
    var idx3  = 0
    var idx4  = 0
    var idx5  = 0
    mh = n >> 3
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) + a(idx2)
    x0i = -a(offa + 1) - a(idx2 + 1)
    x1r = a(offa) - a(idx2)
    x1i = -a(offa + 1) + a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i - x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i + x2i
    a(idx2) = x1r + x3i
    a(idx2 + 1) = x1i + x3r
    a(idx3) = x1r - x3i
    a(idx3 + 1) = x1i - x3r
    wn4r = w(startw + 1)
    csc1 = w(startw + 2)
    csc3 = w(startw + 3)
    wd1r = 1
    wd1i = 0
    wd3r = 1
    wd3i = 0
    k = 0
    var j = 2
    while ( {
      j < mh - 2
    }) {
      k += 4
      idx4 = startw + k
      wk1r = csc1 * (wd1r + w(idx4))
      wk1i = csc1 * (wd1i + w(idx4 + 1))
      wk3r = csc3 * (wd3r + w(idx4 + 2))
      wk3i = csc3 * (wd3i + w(idx4 + 3))
      wd1r = w(idx4)
      wd1i = w(idx4 + 1)
      wd3r = w(idx4 + 2)
      wd3i = w(idx4 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx5 = offa + j
      x0r = a(idx5) + a(idx2)
      x0i = -a(idx5 + 1) - a(idx2 + 1)
      x1r = a(idx5) - a(offa + j2)
      x1i = -a(idx5 + 1) + a(idx2 + 1)
      y0r = a(idx5 + 2) + a(idx2 + 2)
      y0i = -a(idx5 + 3) - a(idx2 + 3)
      y1r = a(idx5 + 2) - a(idx2 + 2)
      y1i = -a(idx5 + 3) + a(idx2 + 3)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 + 2) + a(idx3 + 2)
      y2i = a(idx1 + 3) + a(idx3 + 3)
      y3r = a(idx1 + 2) - a(idx3 + 2)
      y3i = a(idx1 + 3) - a(idx3 + 3)
      a(idx5) = x0r + x2r
      a(idx5 + 1) = x0i - x2i
      a(idx5 + 2) = y0r + y2r
      a(idx5 + 3) = y0i - y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i + x2i
      a(idx1 + 2) = y0r - y2r
      a(idx1 + 3) = y0i + y2i
      x0r = x1r + x3i
      x0i = x1i + x3r
      a(idx2) = wk1r * x0r - wk1i * x0i
      a(idx2 + 1) = wk1r * x0i + wk1i * x0r
      x0r = y1r + y3i
      x0i = y1i + y3r
      a(idx2 + 2) = wd1r * x0r - wd1i * x0i
      a(idx2 + 3) = wd1r * x0i + wd1i * x0r
      x0r = x1r - x3i
      x0i = x1i - x3r
      a(idx3) = wk3r * x0r + wk3i * x0i
      a(idx3 + 1) = wk3r * x0i - wk3i * x0r
      x0r = y1r - y3i
      x0i = y1i - y3r
      a(idx3 + 2) = wd3r * x0r + wd3i * x0i
      a(idx3 + 3) = wd3r * x0i - wd3i * x0r
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) + a(idx2)
      x0i = -a(idx0 + 1) - a(idx2 + 1)
      x1r = a(idx0) - a(idx2)
      x1i = -a(idx0 + 1) + a(idx2 + 1)
      y0r = a(idx0 - 2) + a(idx2 - 2)
      y0i = -a(idx0 - 1) - a(idx2 - 1)
      y1r = a(idx0 - 2) - a(idx2 - 2)
      y1i = -a(idx0 - 1) + a(idx2 - 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      y2r = a(idx1 - 2) + a(idx3 - 2)
      y2i = a(idx1 - 1) + a(idx3 - 1)
      y3r = a(idx1 - 2) - a(idx3 - 2)
      y3i = a(idx1 - 1) - a(idx3 - 1)
      a(idx0) = x0r + x2r
      a(idx0 + 1) = x0i - x2i
      a(idx0 - 2) = y0r + y2r
      a(idx0 - 1) = y0i - y2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i + x2i
      a(idx1 - 2) = y0r - y2r
      a(idx1 - 1) = y0i + y2i
      x0r = x1r + x3i
      x0i = x1i + x3r
      a(idx2) = wk1i * x0r - wk1r * x0i
      a(idx2 + 1) = wk1i * x0i + wk1r * x0r
      x0r = y1r + y3i
      x0i = y1i + y3r
      a(idx2 - 2) = wd1i * x0r - wd1r * x0i
      a(idx2 - 1) = wd1i * x0i + wd1r * x0r
      x0r = x1r - x3i
      x0i = x1i - x3r
      a(idx3) = wk3i * x0r + wk3r * x0i
      a(idx3 + 1) = wk3i * x0i - wk3r * x0r
      x0r = y1r - y3i
      x0i = y1i - y3r
      a(idx3 - 2) = wd3i * x0r + wd3r * x0i
      a(idx3 - 1) = wd3i * x0i - wd3r * x0r

      j += 4
    }
    wk1r = csc1 * (wd1r + wn4r)
    wk1i = csc1 * (wd1i + wn4r)
    wk3r = csc3 * (wd3r - wn4r)
    wk3i = csc3 * (wd3i - wn4r)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0 - 2) + a(idx2 - 2)
    x0i = -a(idx0 - 1) - a(idx2 - 1)
    x1r = a(idx0 - 2) - a(idx2 - 2)
    x1i = -a(idx0 - 1) + a(idx2 - 1)
    x2r = a(idx1 - 2) + a(idx3 - 2)
    x2i = a(idx1 - 1) + a(idx3 - 1)
    x3r = a(idx1 - 2) - a(idx3 - 2)
    x3i = a(idx1 - 1) - a(idx3 - 1)
    a(idx0 - 2) = x0r + x2r
    a(idx0 - 1) = x0i - x2i
    a(idx1 - 2) = x0r - x2r
    a(idx1 - 1) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(idx2 - 2) = wk1r * x0r - wk1i * x0i
    a(idx2 - 1) = wk1r * x0i + wk1i * x0r
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(idx3 - 2) = wk3r * x0r + wk3i * x0i
    a(idx3 - 1) = wk3r * x0i - wk3i * x0r
    x0r = a(idx0) + a(idx2)
    x0i = -a(idx0 + 1) - a(idx2 + 1)
    x1r = a(idx0) - a(idx2)
    x1i = -a(idx0 + 1) + a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(idx0) = x0r + x2r
    a(idx0 + 1) = x0i - x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(idx2) = wn4r * (x0r - x0i)
    a(idx2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(idx3) = -wn4r * (x0r + x0i)
    a(idx3 + 1) = -wn4r * (x0i - x0r)
    x0r = a(idx0 + 2) + a(idx2 + 2)
    x0i = -a(idx0 + 3) - a(idx2 + 3)
    x1r = a(idx0 + 2) - a(idx2 + 2)
    x1i = -a(idx0 + 3) + a(idx2 + 3)
    x2r = a(idx1 + 2) + a(idx3 + 2)
    x2i = a(idx1 + 3) + a(idx3 + 3)
    x3r = a(idx1 + 2) - a(idx3 + 2)
    x3i = a(idx1 + 3) - a(idx3 + 3)
    a(idx0 + 2) = x0r + x2r
    a(idx0 + 3) = x0i - x2i
    a(idx1 + 2) = x0r - x2r
    a(idx1 + 3) = x0i + x2i
    x0r = x1r + x3i
    x0i = x1i + x3r
    a(idx2 + 2) = wk1i * x0r - wk1r * x0i
    a(idx2 + 3) = wk1i * x0i + wk1r * x0r
    x0r = x1r - x3i
    x0i = x1i - x3r
    a(idx3 + 2) = wk3i * x0r + wk3r * x0i
    a(idx3 + 3) = wk3i * x0i - wk3r * x0r
  }

  def cftrec4_th(n: Int, a: Array[Float], offa: Int, nw: Int, w: Array[Float]): Unit = {
    var i = 0
    var idiv4 = 0
    var m = 0
    var nthreads = 0
    var idx = 0
    nthreads = 2
    idiv4 = 0
    m = n >> 1
    if (n >= CommonUtils.threadsBeginN_1D_FFT_4Threads) {
      nthreads = 4
      idiv4 = 1
      m >>= 1
    }
    val futures = new Array[Future[_]](nthreads)
    val mf = m
    i = 0
    while (i < nthreads) {
      val firstIdx = offa + i * m
      futures(idx) = if (i != idiv4) Future {
        var isplt = 0
        var j     = 0
        var k     = 0
        var m     = 0
        val idx1  = firstIdx + mf
        m = n
        while (m > 512) {
          m >>= 2
          cftmdl1(m, a, idx1 - m, w, nw - (m >> 1))
        }
        cftleaf(m, 1, a, idx1 - m, nw, w)
        k = 0
        val idx2 = firstIdx - m
        j = mf - m
        while (j > 0) {
          k += 1
          isplt = cfttree(m, j, k, a, firstIdx, nw, w)
          cftleaf(m, isplt, a, idx2 + j, nw, w)

          j -= m
        }
      }
      else Future {
        var isplt = 0
        var j     = 0
        var k     = 0
        var m     = 0
        val idx1  = firstIdx + mf
        k = 1
        m = n
        while (m > 512) {
          m >>= 2
          k <<= 2
          cftmdl2(m, a, idx1 - m, w, nw - m)
        }
        cftleaf(m, 0, a, idx1 - m, nw, w)
        k >>= 1
        val idx2 = firstIdx - m
        j = mf - m
        while (j > 0) {
          k += 1
          isplt = cfttree(m, j, k, a, firstIdx, nw, w)
          cftleaf(m, isplt, a, idx2 + j, nw, w)

          j -= m
        }
      }
      idx += 1

      i += 1
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  def cftrec4(n: Int, a: Array[Float], offa: Int, nw: Int, w: Array[Float]): Unit = {
    var isplt = 0
    var j = 0
    var k = 0
    var m = 0
    m = n
    val idx1 = offa + n
    while ( {
      m > 512
    }) {
      m >>= 2
      cftmdl1(m, a, idx1 - m, w, nw - (m >> 1))
    }
    cftleaf(m, 1, a, idx1 - m, nw, w)
    k = 0
    val idx2 = offa - m
    j = n - m
    while ( {
      j > 0
    }) {
      k += 1
      isplt = cfttree(m, j, k, a, offa, nw, w)
      cftleaf(m, isplt, a, idx2 + j, nw, w)

      j -= m
    }
  }

  def cfttree(n: Int, j: Int, k: Int, a: Array[Float], offa: Int, nw: Int, w: Array[Float]): Int = {
    var i = 0
    var isplt = 0
    var m = 0
    val idx1 = offa - n
    if ((k & 3) != 0) {
      isplt = k & 1
      if (isplt != 0) cftmdl1(n, a, idx1 + j, w, nw - (n >> 1))
      else cftmdl2(n, a, idx1 + j, w, nw - n)
    }
    else {
      m = n
      i = k
      while ( {
        (i & 3) == 0
      }) {
        m <<= 2

        i >>= 2
      }
      isplt = i & 1
      val idx2 = offa + j
      if (isplt != 0) while ( {
        m > 128
      }) {
        cftmdl1(m, a, idx2 - m, w, nw - (m >> 1))
        m >>= 2
      }
      else while ( {
        m > 128
      }) {
        cftmdl2(m, a, idx2 - m, w, nw - m)
        m >>= 2
      }
    }
    isplt
  }

  def cftleaf(n: Int, isplt: Int, a: Array[Float], offa: Int, nw: Int, w: Array[Float]): Unit = {
    if (n == 512) {
      cftmdl1(128, a, offa, w, nw - 64)
      cftf161(a, offa, w, nw - 8)
      cftf162(a, offa + 32, w, nw - 32)
      cftf161(a, offa + 64, w, nw - 8)
      cftf161(a, offa + 96, w, nw - 8)
      cftmdl2(128, a, offa + 128, w, nw - 128)
      cftf161(a, offa + 128, w, nw - 8)
      cftf162(a, offa + 160, w, nw - 32)
      cftf161(a, offa + 192, w, nw - 8)
      cftf162(a, offa + 224, w, nw - 32)
      cftmdl1(128, a, offa + 256, w, nw - 64)
      cftf161(a, offa + 256, w, nw - 8)
      cftf162(a, offa + 288, w, nw - 32)
      cftf161(a, offa + 320, w, nw - 8)
      cftf161(a, offa + 352, w, nw - 8)
      if (isplt != 0) {
        cftmdl1(128, a, offa + 384, w, nw - 64)
        cftf161(a, offa + 480, w, nw - 8)
      }
      else {
        cftmdl2(128, a, offa + 384, w, nw - 128)
        cftf162(a, offa + 480, w, nw - 32)
      }
      cftf161(a, offa + 384, w, nw - 8)
      cftf162(a, offa + 416, w, nw - 32)
      cftf161(a, offa + 448, w, nw - 8)
    }
    else {
      cftmdl1(64, a, offa, w, nw - 32)
      cftf081(a, offa, w, nw - 8)
      cftf082(a, offa + 16, w, nw - 8)
      cftf081(a, offa + 32, w, nw - 8)
      cftf081(a, offa + 48, w, nw - 8)
      cftmdl2(64, a, offa + 64, w, nw - 64)
      cftf081(a, offa + 64, w, nw - 8)
      cftf082(a, offa + 80, w, nw - 8)
      cftf081(a, offa + 96, w, nw - 8)
      cftf082(a, offa + 112, w, nw - 8)
      cftmdl1(64, a, offa + 128, w, nw - 32)
      cftf081(a, offa + 128, w, nw - 8)
      cftf082(a, offa + 144, w, nw - 8)
      cftf081(a, offa + 160, w, nw - 8)
      cftf081(a, offa + 176, w, nw - 8)
      if (isplt != 0) {
        cftmdl1(64, a, offa + 192, w, nw - 32)
        cftf081(a, offa + 240, w, nw - 8)
      }
      else {
        cftmdl2(64, a, offa + 192, w, nw - 64)
        cftf082(a, offa + 240, w, nw - 8)
      }
      cftf081(a, offa + 192, w, nw - 8)
      cftf082(a, offa + 208, w, nw - 8)
      cftf081(a, offa + 224, w, nw - 8)
    }
  }

  def cftmdl1(n: Int, a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var j0    = 0
    var j1    = 0
    var j2    = 0
    var j3    = 0
    var k     = 0
    var m     = 0
    var mh    = 0
    var wn4r  = 0.0f
    var wk1r  = 0.0f
    var wk1i  = 0.0f
    var wk3r  = 0.0f
    var wk3i  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var idx0  = 0
    var idx1  = 0
    var idx2  = 0
    var idx3  = 0
    var idx4  = 0
    var idx5  = 0
    mh = n >> 3
    m = 2 * mh
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) + a(idx2)
    x0i = a(offa + 1) + a(idx2 + 1)
    x1r = a(offa) - a(idx2)
    x1i = a(offa + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    a(idx2) = x1r - x3i
    a(idx2 + 1) = x1i + x3r
    a(idx3) = x1r + x3i
    a(idx3 + 1) = x1i - x3r
    wn4r = w(startw + 1)
    k = 0
    var j = 2
    while ( {
      j < mh
    }) {
      k += 4
      idx4 = startw + k
      wk1r = w(idx4)
      wk1i = w(idx4 + 1)
      wk3r = w(idx4 + 2)
      wk3i = w(idx4 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx5 = offa + j
      x0r = a(idx5) + a(idx2)
      x0i = a(idx5 + 1) + a(idx2 + 1)
      x1r = a(idx5) - a(idx2)
      x1i = a(idx5 + 1) - a(idx2 + 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      a(idx5) = x0r + x2r
      a(idx5 + 1) = x0i + x2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1r * x0r - wk1i * x0i
      a(idx2 + 1) = wk1r * x0i + wk1i * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3r * x0r + wk3i * x0i
      a(idx3 + 1) = wk3r * x0i - wk3i * x0r
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) + a(idx2)
      x0i = a(idx0 + 1) + a(idx2 + 1)
      x1r = a(idx0) - a(idx2)
      x1i = a(idx0 + 1) - a(idx2 + 1)
      x2r = a(idx1) + a(idx3)
      x2i = a(idx1 + 1) + a(idx3 + 1)
      x3r = a(idx1) - a(idx3)
      x3i = a(idx1 + 1) - a(idx3 + 1)
      a(idx0) = x0r + x2r
      a(idx0 + 1) = x0i + x2i
      a(idx1) = x0r - x2r
      a(idx1 + 1) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(idx2) = wk1i * x0r - wk1r * x0i
      a(idx2 + 1) = wk1i * x0i + wk1r * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(idx3) = wk3i * x0r + wk3r * x0i
      a(idx3 + 1) = wk3i * x0i - wk3r * x0r

      j += 2
    }
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0) + a(idx2)
    x0i = a(idx0 + 1) + a(idx2 + 1)
    x1r = a(idx0) - a(idx2)
    x1i = a(idx0 + 1) - a(idx2 + 1)
    x2r = a(idx1) + a(idx3)
    x2i = a(idx1 + 1) + a(idx3 + 1)
    x3r = a(idx1) - a(idx3)
    x3i = a(idx1 + 1) - a(idx3 + 1)
    a(idx0) = x0r + x2r
    a(idx0 + 1) = x0i + x2i
    a(idx1) = x0r - x2r
    a(idx1 + 1) = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    a(idx2) = wn4r * (x0r - x0i)
    a(idx2 + 1) = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    a(idx3) = -wn4r * (x0r + x0i)
    a(idx3 + 1) = -wn4r * (x0i - x0r)
  }

  def cftmdl2(n: Int, a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var j0    = 0
    var j1    = 0
    var j2    = 0
    var j3    = 0
    var k     = 0
    var kr    = 0
    var m     = 0
    var mh    = 0
    var wn4r  = 0.0f
    var wk1r  = 0.0f
    var wk1i  = 0.0f
    var wk3r  = 0.0f
    var wk3i  = 0.0f
    var wd1r  = 0.0f
    var wd1i  = 0.0f
    var wd3r  = 0.0f
    var wd3i  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var y0r   = 0.0f
    var y0i   = 0.0f
    var y2r   = 0.0f
    var y2i   = 0.0f
    var idx0  = 0
    var idx1  = 0
    var idx2  = 0
    var idx3  = 0
    var idx4  = 0
    var idx5  = 0
    var idx6  = 0
    mh = n >> 3
    m = 2 * mh
    wn4r = w(startw + 1)
    j1 = m
    j2 = j1 + m
    j3 = j2 + m
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(offa) - a(idx2 + 1)
    x0i = a(offa + 1) + a(idx2)
    x1r = a(offa) + a(idx2 + 1)
    x1i = a(offa + 1) - a(idx2)
    x2r = a(idx1) - a(idx3 + 1)
    x2i = a(idx1 + 1) + a(idx3)
    x3r = a(idx1) + a(idx3 + 1)
    x3i = a(idx1 + 1) - a(idx3)
    y0r = wn4r * (x2r - x2i)
    y0i = wn4r * (x2i + x2r)
    a(offa) = x0r + y0r
    a(offa + 1) = x0i + y0i
    a(idx1) = x0r - y0r
    a(idx1 + 1) = x0i - y0i
    y0r = wn4r * (x3r - x3i)
    y0i = wn4r * (x3i + x3r)
    a(idx2) = x1r - y0i
    a(idx2 + 1) = x1i + y0r
    a(idx3) = x1r + y0i
    a(idx3 + 1) = x1i - y0r
    k = 0
    kr = 2 * m
    var j = 2
    while ( {
      j < mh
    }) {
      k += 4
      idx4 = startw + k
      wk1r = w(idx4)
      wk1i = w(idx4 + 1)
      wk3r = w(idx4 + 2)
      wk3i = w(idx4 + 3)
      kr -= 4
      idx5 = startw + kr
      wd1i = w(idx5)
      wd1r = w(idx5 + 1)
      wd3i = w(idx5 + 2)
      wd3r = w(idx5 + 3)
      j1 = j + m
      j2 = j1 + m
      j3 = j2 + m
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      idx6 = offa + j
      x0r = a(idx6) - a(idx2 + 1)
      x0i = a(idx6 + 1) + a(idx2)
      x1r = a(idx6) + a(idx2 + 1)
      x1i = a(idx6 + 1) - a(idx2)
      x2r = a(idx1) - a(idx3 + 1)
      x2i = a(idx1 + 1) + a(idx3)
      x3r = a(idx1) + a(idx3 + 1)
      x3i = a(idx1 + 1) - a(idx3)
      y0r = wk1r * x0r - wk1i * x0i
      y0i = wk1r * x0i + wk1i * x0r
      y2r = wd1r * x2r - wd1i * x2i
      y2i = wd1r * x2i + wd1i * x2r
      a(idx6) = y0r + y2r
      a(idx6 + 1) = y0i + y2i
      a(idx1) = y0r - y2r
      a(idx1 + 1) = y0i - y2i
      y0r = wk3r * x1r + wk3i * x1i
      y0i = wk3r * x1i - wk3i * x1r
      y2r = wd3r * x3r + wd3i * x3i
      y2i = wd3r * x3i - wd3i * x3r
      a(idx2) = y0r + y2r
      a(idx2 + 1) = y0i + y2i
      a(idx3) = y0r - y2r
      a(idx3 + 1) = y0i - y2i
      j0 = m - j
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      idx0 = offa + j0
      idx1 = offa + j1
      idx2 = offa + j2
      idx3 = offa + j3
      x0r = a(idx0) - a(idx2 + 1)
      x0i = a(idx0 + 1) + a(idx2)
      x1r = a(idx0) + a(idx2 + 1)
      x1i = a(idx0 + 1) - a(idx2)
      x2r = a(idx1) - a(idx3 + 1)
      x2i = a(idx1 + 1) + a(idx3)
      x3r = a(idx1) + a(idx3 + 1)
      x3i = a(idx1 + 1) - a(idx3)
      y0r = wd1i * x0r - wd1r * x0i
      y0i = wd1i * x0i + wd1r * x0r
      y2r = wk1i * x2r - wk1r * x2i
      y2i = wk1i * x2i + wk1r * x2r
      a(idx0) = y0r + y2r
      a(idx0 + 1) = y0i + y2i
      a(idx1) = y0r - y2r
      a(idx1 + 1) = y0i - y2i
      y0r = wd3i * x1r + wd3r * x1i
      y0i = wd3i * x1i - wd3r * x1r
      y2r = wk3i * x3r + wk3r * x3i
      y2i = wk3i * x3i - wk3r * x3r
      a(idx2) = y0r + y2r
      a(idx2 + 1) = y0i + y2i
      a(idx3) = y0r - y2r
      a(idx3 + 1) = y0i - y2i

      j += 2
    }
    wk1r = w(startw + m)
    wk1i = w(startw + m + 1)
    j0 = mh
    j1 = j0 + m
    j2 = j1 + m
    j3 = j2 + m
    idx0 = offa + j0
    idx1 = offa + j1
    idx2 = offa + j2
    idx3 = offa + j3
    x0r = a(idx0) - a(idx2 + 1)
    x0i = a(idx0 + 1) + a(idx2)
    x1r = a(idx0) + a(idx2 + 1)
    x1i = a(idx0 + 1) - a(idx2)
    x2r = a(idx1) - a(idx3 + 1)
    x2i = a(idx1 + 1) + a(idx3)
    x3r = a(idx1) + a(idx3 + 1)
    x3i = a(idx1 + 1) - a(idx3)
    y0r = wk1r * x0r - wk1i * x0i
    y0i = wk1r * x0i + wk1i * x0r
    y2r = wk1i * x2r - wk1r * x2i
    y2i = wk1i * x2i + wk1r * x2r
    a(idx0) = y0r + y2r
    a(idx0 + 1) = y0i + y2i
    a(idx1) = y0r - y2r
    a(idx1 + 1) = y0i - y2i
    y0r = wk1i * x1r - wk1r * x1i
    y0i = wk1i * x1i + wk1r * x1r
    y2r = wk1r * x3r - wk1i * x3i
    y2i = wk1r * x3i + wk1i * x3r
    a(idx2) = y0r - y2r
    a(idx2 + 1) = y0i - y2i
    a(idx3) = y0r + y2r
    a(idx3 + 1) = y0i + y2i
  }

  def cftfx41(n: Int, a: Array[Float], offa: Int, nw: Int, w: Array[Float]): Unit = {
    if (n == 128) {
      cftf161(a, offa, w, nw - 8)
      cftf162(a, offa + 32, w, nw - 32)
      cftf161(a, offa + 64, w, nw - 8)
      cftf161(a, offa + 96, w, nw - 8)
    }
    else {
      cftf081(a, offa, w, nw - 8)
      cftf082(a, offa + 16, w, nw - 8)
      cftf081(a, offa + 32, w, nw - 8)
      cftf081(a, offa + 48, w, nw - 8)
    }
  }

  def cftf161(a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var wn4r  = 0.0f
    var wk1r  = 0.0f
    var wk1i  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var y0r   = 0.0f
    var y0i   = 0.0f
    var y1r   = 0.0f
    var y1i   = 0.0f
    var y2r   = 0.0f
    var y2i   = 0.0f
    var y3r   = 0.0f
    var y3i   = 0.0f
    var y4r   = 0.0f
    var y4i   = 0.0f
    var y5r   = 0.0f
    var y5i   = 0.0f
    var y6r   = 0.0f
    var y6i   = 0.0f
    var y7r   = 0.0f
    var y7i   = 0.0f
    var y8r   = 0.0f
    var y8i   = 0.0f
    var y9r   = 0.0f
    var y9i   = 0.0f
    var y10r  = 0.0f
    var y10i  = 0.0f
    var y11r  = 0.0f
    var y11i  = 0.0f
    var y12r  = 0.0f
    var y12i  = 0.0f
    var y13r  = 0.0f
    var y13i  = 0.0f
    var y14r  = 0.0f
    var y14i  = 0.0f
    var y15r  = 0.0f
    var y15i  = 0.0f
    wn4r = w(startw + 1)
    wk1r = w(startw + 2)
    wk1i = w(startw + 3)
    x0r = a(offa) + a(offa + 16)
    x0i = a(offa + 1) + a(offa + 17)
    x1r = a(offa) - a(offa + 16)
    x1i = a(offa + 1) - a(offa + 17)
    x2r = a(offa + 8) + a(offa + 24)
    x2i = a(offa + 9) + a(offa + 25)
    x3r = a(offa + 8) - a(offa + 24)
    x3i = a(offa + 9) - a(offa + 25)
    y0r = x0r + x2r
    y0i = x0i + x2i
    y4r = x0r - x2r
    y4i = x0i - x2i
    y8r = x1r - x3i
    y8i = x1i + x3r
    y12r = x1r + x3i
    y12i = x1i - x3r
    x0r = a(offa + 2) + a(offa + 18)
    x0i = a(offa + 3) + a(offa + 19)
    x1r = a(offa + 2) - a(offa + 18)
    x1i = a(offa + 3) - a(offa + 19)
    x2r = a(offa + 10) + a(offa + 26)
    x2i = a(offa + 11) + a(offa + 27)
    x3r = a(offa + 10) - a(offa + 26)
    x3i = a(offa + 11) - a(offa + 27)
    y1r = x0r + x2r
    y1i = x0i + x2i
    y5r = x0r - x2r
    y5i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y9r = wk1r * x0r - wk1i * x0i
    y9i = wk1r * x0i + wk1i * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    y13r = wk1i * x0r - wk1r * x0i
    y13i = wk1i * x0i + wk1r * x0r
    x0r = a(offa + 4) + a(offa + 20)
    x0i = a(offa + 5) + a(offa + 21)
    x1r = a(offa + 4) - a(offa + 20)
    x1i = a(offa + 5) - a(offa + 21)
    x2r = a(offa + 12) + a(offa + 28)
    x2i = a(offa + 13) + a(offa + 29)
    x3r = a(offa + 12) - a(offa + 28)
    x3i = a(offa + 13) - a(offa + 29)
    y2r = x0r + x2r
    y2i = x0i + x2i
    y6r = x0r - x2r
    y6i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y10r = wn4r * (x0r - x0i)
    y10i = wn4r * (x0i + x0r)
    x0r = x1r + x3i
    x0i = x1i - x3r
    y14r = wn4r * (x0r + x0i)
    y14i = wn4r * (x0i - x0r)
    x0r = a(offa + 6) + a(offa + 22)
    x0i = a(offa + 7) + a(offa + 23)
    x1r = a(offa + 6) - a(offa + 22)
    x1i = a(offa + 7) - a(offa + 23)
    x2r = a(offa + 14) + a(offa + 30)
    x2i = a(offa + 15) + a(offa + 31)
    x3r = a(offa + 14) - a(offa + 30)
    x3i = a(offa + 15) - a(offa + 31)
    y3r = x0r + x2r
    y3i = x0i + x2i
    y7r = x0r - x2r
    y7i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    y11r = wk1i * x0r - wk1r * x0i
    y11i = wk1i * x0i + wk1r * x0r
    x0r = x1r + x3i
    x0i = x1i - x3r
    y15r = wk1r * x0r - wk1i * x0i
    y15i = wk1r * x0i + wk1i * x0r
    x0r = y12r - y14r
    x0i = y12i - y14i
    x1r = y12r + y14r
    x1i = y12i + y14i
    x2r = y13r - y15r
    x2i = y13i - y15i
    x3r = y13r + y15r
    x3i = y13i + y15i
    a(offa + 24) = x0r + x2r
    a(offa + 25) = x0i + x2i
    a(offa + 26) = x0r - x2r
    a(offa + 27) = x0i - x2i
    a(offa + 28) = x1r - x3i
    a(offa + 29) = x1i + x3r
    a(offa + 30) = x1r + x3i
    a(offa + 31) = x1i - x3r
    x0r = y8r + y10r
    x0i = y8i + y10i
    x1r = y8r - y10r
    x1i = y8i - y10i
    x2r = y9r + y11r
    x2i = y9i + y11i
    x3r = y9r - y11r
    x3i = y9i - y11i
    a(offa + 16) = x0r + x2r
    a(offa + 17) = x0i + x2i
    a(offa + 18) = x0r - x2r
    a(offa + 19) = x0i - x2i
    a(offa + 20) = x1r - x3i
    a(offa + 21) = x1i + x3r
    a(offa + 22) = x1r + x3i
    a(offa + 23) = x1i - x3r
    x0r = y5r - y7i
    x0i = y5i + y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    x0r = y5r + y7i
    x0i = y5i - y7r
    x3r = wn4r * (x0r - x0i)
    x3i = wn4r * (x0i + x0r)
    x0r = y4r - y6i
    x0i = y4i + y6r
    x1r = y4r + y6i
    x1i = y4i - y6r
    a(offa + 8) = x0r + x2r
    a(offa + 9) = x0i + x2i
    a(offa + 10) = x0r - x2r
    a(offa + 11) = x0i - x2i
    a(offa + 12) = x1r - x3i
    a(offa + 13) = x1i + x3r
    a(offa + 14) = x1r + x3i
    a(offa + 15) = x1i - x3r
    x0r = y0r + y2r
    x0i = y0i + y2i
    x1r = y0r - y2r
    x1i = y0i - y2i
    x2r = y1r + y3r
    x2i = y1i + y3i
    x3r = y1r - y3r
    x3i = y1i - y3i
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(offa + 2) = x0r - x2r
    a(offa + 3) = x0i - x2i
    a(offa + 4) = x1r - x3i
    a(offa + 5) = x1i + x3r
    a(offa + 6) = x1r + x3i
    a(offa + 7) = x1i - x3r
  }

  def cftf162(a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var wn4r  = 0.0f
    var wk1r  = 0.0f
    var wk1i  = 0.0f
    var wk2r  = 0.0f
    var wk2i  = 0.0f
    var wk3r  = 0.0f
    var wk3i  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var y0r   = 0.0f
    var y0i   = 0.0f
    var y1r   = 0.0f
    var y1i   = 0.0f
    var y2r   = 0.0f
    var y2i   = 0.0f
    var y3r   = 0.0f
    var y3i   = 0.0f
    var y4r   = 0.0f
    var y4i   = 0.0f
    var y5r   = 0.0f
    var y5i   = 0.0f
    var y6r   = 0.0f
    var y6i   = 0.0f
    var y7r   = 0.0f
    var y7i   = 0.0f
    var y8r   = 0.0f
    var y8i   = 0.0f
    var y9r   = 0.0f
    var y9i   = 0.0f
    var y10r  = 0.0f
    var y10i  = 0.0f
    var y11r  = 0.0f
    var y11i  = 0.0f
    var y12r  = 0.0f
    var y12i  = 0.0f
    var y13r  = 0.0f
    var y13i  = 0.0f
    var y14r  = 0.0f
    var y14i  = 0.0f
    var y15r  = 0.0f
    var y15i  = 0.0f
    wn4r = w(startw + 1)
    wk1r = w(startw + 4)
    wk1i = w(startw + 5)
    wk3r = w(startw + 6)
    wk3i = -w(startw + 7)
    wk2r = w(startw + 8)
    wk2i = w(startw + 9)
    x1r = a(offa) - a(offa + 17)
    x1i = a(offa + 1) + a(offa + 16)
    x0r = a(offa + 8) - a(offa + 25)
    x0i = a(offa + 9) + a(offa + 24)
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    y0r = x1r + x2r
    y0i = x1i + x2i
    y4r = x1r - x2r
    y4i = x1i - x2i
    x1r = a(offa) + a(offa + 17)
    x1i = a(offa + 1) - a(offa + 16)
    x0r = a(offa + 8) + a(offa + 25)
    x0i = a(offa + 9) - a(offa + 24)
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    y8r = x1r - x2i
    y8i = x1i + x2r
    y12r = x1r + x2i
    y12i = x1i - x2r
    x0r = a(offa + 2) - a(offa + 19)
    x0i = a(offa + 3) + a(offa + 18)
    x1r = wk1r * x0r - wk1i * x0i
    x1i = wk1r * x0i + wk1i * x0r
    x0r = a(offa + 10) - a(offa + 27)
    x0i = a(offa + 11) + a(offa + 26)
    x2r = wk3i * x0r - wk3r * x0i
    x2i = wk3i * x0i + wk3r * x0r
    y1r = x1r + x2r
    y1i = x1i + x2i
    y5r = x1r - x2r
    y5i = x1i - x2i
    x0r = a(offa + 2) + a(offa + 19)
    x0i = a(offa + 3) - a(offa + 18)
    x1r = wk3r * x0r - wk3i * x0i
    x1i = wk3r * x0i + wk3i * x0r
    x0r = a(offa + 10) + a(offa + 27)
    x0i = a(offa + 11) - a(offa + 26)
    x2r = wk1r * x0r + wk1i * x0i
    x2i = wk1r * x0i - wk1i * x0r
    y9r = x1r - x2r
    y9i = x1i - x2i
    y13r = x1r + x2r
    y13i = x1i + x2i
    x0r = a(offa + 4) - a(offa + 21)
    x0i = a(offa + 5) + a(offa + 20)
    x1r = wk2r * x0r - wk2i * x0i
    x1i = wk2r * x0i + wk2i * x0r
    x0r = a(offa + 12) - a(offa + 29)
    x0i = a(offa + 13) + a(offa + 28)
    x2r = wk2i * x0r - wk2r * x0i
    x2i = wk2i * x0i + wk2r * x0r
    y2r = x1r + x2r
    y2i = x1i + x2i
    y6r = x1r - x2r
    y6i = x1i - x2i
    x0r = a(offa + 4) + a(offa + 21)
    x0i = a(offa + 5) - a(offa + 20)
    x1r = wk2i * x0r - wk2r * x0i
    x1i = wk2i * x0i + wk2r * x0r
    x0r = a(offa + 12) + a(offa + 29)
    x0i = a(offa + 13) - a(offa + 28)
    x2r = wk2r * x0r - wk2i * x0i
    x2i = wk2r * x0i + wk2i * x0r
    y10r = x1r - x2r
    y10i = x1i - x2i
    y14r = x1r + x2r
    y14i = x1i + x2i
    x0r = a(offa + 6) - a(offa + 23)
    x0i = a(offa + 7) + a(offa + 22)
    x1r = wk3r * x0r - wk3i * x0i
    x1i = wk3r * x0i + wk3i * x0r
    x0r = a(offa + 14) - a(offa + 31)
    x0i = a(offa + 15) + a(offa + 30)
    x2r = wk1i * x0r - wk1r * x0i
    x2i = wk1i * x0i + wk1r * x0r
    y3r = x1r + x2r
    y3i = x1i + x2i
    y7r = x1r - x2r
    y7i = x1i - x2i
    x0r = a(offa + 6) + a(offa + 23)
    x0i = a(offa + 7) - a(offa + 22)
    x1r = wk1i * x0r + wk1r * x0i
    x1i = wk1i * x0i - wk1r * x0r
    x0r = a(offa + 14) + a(offa + 31)
    x0i = a(offa + 15) - a(offa + 30)
    x2r = wk3i * x0r - wk3r * x0i
    x2i = wk3i * x0i + wk3r * x0r
    y11r = x1r + x2r
    y11i = x1i + x2i
    y15r = x1r - x2r
    y15i = x1i - x2i
    x1r = y0r + y2r
    x1i = y0i + y2i
    x2r = y1r + y3r
    x2i = y1i + y3i
    a(offa) = x1r + x2r
    a(offa + 1) = x1i + x2i
    a(offa + 2) = x1r - x2r
    a(offa + 3) = x1i - x2i
    x1r = y0r - y2r
    x1i = y0i - y2i
    x2r = y1r - y3r
    x2i = y1i - y3i
    a(offa + 4) = x1r - x2i
    a(offa + 5) = x1i + x2r
    a(offa + 6) = x1r + x2i
    a(offa + 7) = x1i - x2r
    x1r = y4r - y6i
    x1i = y4i + y6r
    x0r = y5r - y7i
    x0i = y5i + y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 8) = x1r + x2r
    a(offa + 9) = x1i + x2i
    a(offa + 10) = x1r - x2r
    a(offa + 11) = x1i - x2i
    x1r = y4r + y6i
    x1i = y4i - y6r
    x0r = y5r + y7i
    x0i = y5i - y7r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 12) = x1r - x2i
    a(offa + 13) = x1i + x2r
    a(offa + 14) = x1r + x2i
    a(offa + 15) = x1i - x2r
    x1r = y8r + y10r
    x1i = y8i + y10i
    x2r = y9r - y11r
    x2i = y9i - y11i
    a(offa + 16) = x1r + x2r
    a(offa + 17) = x1i + x2i
    a(offa + 18) = x1r - x2r
    a(offa + 19) = x1i - x2i
    x1r = y8r - y10r
    x1i = y8i - y10i
    x2r = y9r + y11r
    x2i = y9i + y11i
    a(offa + 20) = x1r - x2i
    a(offa + 21) = x1i + x2r
    a(offa + 22) = x1r + x2i
    a(offa + 23) = x1i - x2r
    x1r = y12r - y14i
    x1i = y12i + y14r
    x0r = y13r + y15i
    x0i = y13i - y15r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 24) = x1r + x2r
    a(offa + 25) = x1i + x2i
    a(offa + 26) = x1r - x2r
    a(offa + 27) = x1i - x2i
    x1r = y12r + y14i
    x1i = y12i - y14r
    x0r = y13r - y15i
    x0i = y13i + y15r
    x2r = wn4r * (x0r - x0i)
    x2i = wn4r * (x0i + x0r)
    a(offa + 28) = x1r - x2i
    a(offa + 29) = x1i + x2r
    a(offa + 30) = x1r + x2i
    a(offa + 31) = x1i - x2r
  }

  def cftf081(a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var wn4r  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var x2r   = 0.0f
    var x2i   = 0.0f
    var x3r   = 0.0f
    var x3i   = 0.0f
    var y0r   = 0.0f
    var y0i   = 0.0f
    var y1r   = 0.0f
    var y1i   = 0.0f
    var y2r   = 0.0f
    var y2i   = 0.0f
    var y3r   = 0.0f
    var y3i   = 0.0f
    var y4r   = 0.0f
    var y4i   = 0.0f
    var y5r   = 0.0f
    var y5i   = 0.0f
    var y6r   = 0.0f
    var y6i   = 0.0f
    var y7r   = 0.0f
    var y7i   = 0.0f
    wn4r = w(startw + 1)
    x0r = a(offa) + a(offa + 8)
    x0i = a(offa + 1) + a(offa + 9)
    x1r = a(offa) - a(offa + 8)
    x1i = a(offa + 1) - a(offa + 9)
    x2r = a(offa + 4) + a(offa + 12)
    x2i = a(offa + 5) + a(offa + 13)
    x3r = a(offa + 4) - a(offa + 12)
    x3i = a(offa + 5) - a(offa + 13)
    y0r = x0r + x2r
    y0i = x0i + x2i
    y2r = x0r - x2r
    y2i = x0i - x2i
    y1r = x1r - x3i
    y1i = x1i + x3r
    y3r = x1r + x3i
    y3i = x1i - x3r
    x0r = a(offa + 2) + a(offa + 10)
    x0i = a(offa + 3) + a(offa + 11)
    x1r = a(offa + 2) - a(offa + 10)
    x1i = a(offa + 3) - a(offa + 11)
    x2r = a(offa + 6) + a(offa + 14)
    x2i = a(offa + 7) + a(offa + 15)
    x3r = a(offa + 6) - a(offa + 14)
    x3i = a(offa + 7) - a(offa + 15)
    y4r = x0r + x2r
    y4i = x0i + x2i
    y6r = x0r - x2r
    y6i = x0i - x2i
    x0r = x1r - x3i
    x0i = x1i + x3r
    x2r = x1r + x3i
    x2i = x1i - x3r
    y5r = wn4r * (x0r - x0i)
    y5i = wn4r * (x0r + x0i)
    y7r = wn4r * (x2r - x2i)
    y7i = wn4r * (x2r + x2i)
    a(offa + 8) = y1r + y5r
    a(offa + 9) = y1i + y5i
    a(offa + 10) = y1r - y5r
    a(offa + 11) = y1i - y5i
    a(offa + 12) = y3r - y7i
    a(offa + 13) = y3i + y7r
    a(offa + 14) = y3r + y7i
    a(offa + 15) = y3i - y7r
    a(offa) = y0r + y4r
    a(offa + 1) = y0i + y4i
    a(offa + 2) = y0r - y4r
    a(offa + 3) = y0i - y4i
    a(offa + 4) = y2r - y6i
    a(offa + 5) = y2i + y6r
    a(offa + 6) = y2r + y6i
    a(offa + 7) = y2i - y6r
  }

  def cftf082(a: Array[Float], offa: Int, w: Array[Float], startw: Int): Unit = {
    var wn4r  = 0.0f
    var wk1r  = 0.0f
    var wk1i  = 0.0f
    var x0r   = 0.0f
    var x0i   = 0.0f
    var x1r   = 0.0f
    var x1i   = 0.0f
    var y0r   = 0.0f
    var y0i   = 0.0f
    var y1r   = 0.0f
    var y1i   = 0.0f
    var y2r   = 0.0f
    var y2i   = 0.0f
    var y3r   = 0.0f
    var y3i   = 0.0f
    var y4r   = 0.0f
    var y4i   = 0.0f
    var y5r   = 0.0f
    var y5i   = 0.0f
    var y6r   = 0.0f
    var y6i   = 0.0f
    var y7r   = 0.0f
    var y7i   = 0.0f
    wn4r = w(startw + 1)
    wk1r = w(startw + 2)
    wk1i = w(startw + 3)
    y0r = a(offa) - a(offa + 9)
    y0i = a(offa + 1) + a(offa + 8)
    y1r = a(offa) + a(offa + 9)
    y1i = a(offa + 1) - a(offa + 8)
    x0r = a(offa + 4) - a(offa + 13)
    x0i = a(offa + 5) + a(offa + 12)
    y2r = wn4r * (x0r - x0i)
    y2i = wn4r * (x0i + x0r)
    x0r = a(offa + 4) + a(offa + 13)
    x0i = a(offa + 5) - a(offa + 12)
    y3r = wn4r * (x0r - x0i)
    y3i = wn4r * (x0i + x0r)
    x0r = a(offa + 2) - a(offa + 11)
    x0i = a(offa + 3) + a(offa + 10)
    y4r = wk1r * x0r - wk1i * x0i
    y4i = wk1r * x0i + wk1i * x0r
    x0r = a(offa + 2) + a(offa + 11)
    x0i = a(offa + 3) - a(offa + 10)
    y5r = wk1i * x0r - wk1r * x0i
    y5i = wk1i * x0i + wk1r * x0r
    x0r = a(offa + 6) - a(offa + 15)
    x0i = a(offa + 7) + a(offa + 14)
    y6r = wk1i * x0r - wk1r * x0i
    y6i = wk1i * x0i + wk1r * x0r
    x0r = a(offa + 6) + a(offa + 15)
    x0i = a(offa + 7) - a(offa + 14)
    y7r = wk1r * x0r - wk1i * x0i
    y7i = wk1r * x0i + wk1i * x0r
    x0r = y0r + y2r
    x0i = y0i + y2i
    x1r = y4r + y6r
    x1i = y4i + y6i
    a(offa) = x0r + x1r
    a(offa + 1) = x0i + x1i
    a(offa + 2) = x0r - x1r
    a(offa + 3) = x0i - x1i
    x0r = y0r - y2r
    x0i = y0i - y2i
    x1r = y4r - y6r
    x1i = y4i - y6i
    a(offa + 4) = x0r - x1i
    a(offa + 5) = x0i + x1r
    a(offa + 6) = x0r + x1i
    a(offa + 7) = x0i - x1r
    x0r = y1r - y3i
    x0i = y1i + y3r
    x1r = y5r - y7r
    x1i = y5i - y7i
    a(offa + 8) = x0r + x1r
    a(offa + 9) = x0i + x1i
    a(offa + 10) = x0r - x1r
    a(offa + 11) = x0i - x1i
    x0r = y1r + y3i
    x0i = y1i - y3r
    x1r = y5r + y7r
    x1i = y5i + y7i
    a(offa + 12) = x0r - x1i
    a(offa + 13) = x0i + x1r
    a(offa + 14) = x0r + x1i
    a(offa + 15) = x0i - x1r
  }

  def cftf040(a: Array[Float], offa: Int): Unit = {
    var x0r = 0.0f
    var x0i = 0.0f
    var x1r = 0.0f
    var x1i = 0.0f
    var x2r = 0.0f
    var x2i = 0.0f
    var x3r = 0.0f
    var x3i = 0.0f
    x0r = a(offa) + a(offa + 4)
    x0i = a(offa + 1) + a(offa + 5)
    x1r = a(offa) - a(offa + 4)
    x1i = a(offa + 1) - a(offa + 5)
    x2r = a(offa + 2) + a(offa + 6)
    x2i = a(offa + 3) + a(offa + 7)
    x3r = a(offa + 2) - a(offa + 6)
    x3i = a(offa + 3) - a(offa + 7)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(offa + 2) = x1r - x3i
    a(offa + 3) = x1i + x3r
    a(offa + 4) = x0r - x2r
    a(offa + 5) = x0i - x2i
    a(offa + 6) = x1r + x3i
    a(offa + 7) = x1i - x3r
  }

  def cftb040(a: Array[Float], offa: Int): Unit = {
    var x0r = 0.0f
    var x0i = 0.0f
    var x1r = 0.0f
    var x1i = 0.0f
    var x2r = 0.0f
    var x2i = 0.0f
    var x3r = 0.0f
    var x3i = 0.0f
    x0r = a(offa) + a(offa + 4)
    x0i = a(offa + 1) + a(offa + 5)
    x1r = a(offa) - a(offa + 4)
    x1i = a(offa + 1) - a(offa + 5)
    x2r = a(offa + 2) + a(offa + 6)
    x2i = a(offa + 3) + a(offa + 7)
    x3r = a(offa + 2) - a(offa + 6)
    x3i = a(offa + 3) - a(offa + 7)
    a(offa) = x0r + x2r
    a(offa + 1) = x0i + x2i
    a(offa + 2) = x1r + x3i
    a(offa + 3) = x1i - x3r
    a(offa + 4) = x0r - x2r
    a(offa + 5) = x0i - x2i
    a(offa + 6) = x1r - x3i
    a(offa + 7) = x1i + x3r
  }

  def cftx020(a: Array[Float], offa: Int): Unit = {
    var x0r = 0.0f
    var x0i = 0.0f
    x0r = a(offa) - a(offa + 2)
    x0i = -a(offa + 1) + a(offa + 3)
    a(offa) += a(offa + 2)
    a(offa + 1) += a(offa + 3)
    a(offa + 2) = x0r
    a(offa + 3) = x0i
  }

  def cftxb020(a: Array[Float], offa: Int): Unit = {
    var x0r = 0.0f
    var x0i = 0.0f
    x0r = a(offa) - a(offa + 2)
    x0i = a(offa + 1) - a(offa + 3)
    a(offa) += a(offa + 2)
    a(offa + 1) += a(offa + 3)
    a(offa + 2) = x0r
    a(offa + 3) = x0i
  }

  def cftxc020(a: Array[Float], offa: Int): Unit = {
    var x0r = 0.0f
    var x0i = 0.0f
    x0r = a(offa) - a(offa + 2)
    x0i = a(offa + 1) + a(offa + 3)
    a(offa) += a(offa + 2)
    a(offa + 1) -= a(offa + 3)
    a(offa + 2) = x0r
    a(offa + 3) = x0i
  }

  def rftfsub(n: Int, a: Array[Float], offa: Int, nc: Int, c: Array[Float], startc: Int): Unit = {
    var k     = 0
    var kk    = 0
    var ks    = 0
    var m     = 0
    var wkr   = 0.0f
    var wki   = 0.0f
    var xr    = 0.0f
    var xi    = 0.0f
    var yr    = 0.0f
    var yi    = 0.0f
    var idx1  = 0
    var idx2  = 0
    m = n >> 1
    ks = 2 * nc / m
    kk = 0
    var j = 2
    while ( {
      j < m
    }) {
      k = n - j
      kk += ks
      wkr = 0.5f - c(startc + nc - kk)
      wki = c(startc + kk)
      idx1 = offa + j
      idx2 = offa + k
      xr = a(idx1) - a(idx2)
      xi = a(idx1 + 1) + a(idx2 + 1)
      yr = wkr * xr - wki * xi
      yi = wkr * xi + wki * xr
      a(idx1) -= yr
      a(idx1 + 1) = yi - a(idx1 + 1)
      a(idx2) += yr
      a(idx2 + 1) = yi - a(idx2 + 1)

      j += 2
    }
    a(offa + m + 1) = -a(offa + m + 1)
  }

  def rftbsub(n: Int, a: Array[Float], offa: Int, nc: Int, c: Array[Float], startc: Int): Unit = {
    var k     = 0
    var kk    = 0
    var ks    = 0
    var m     = 0
    var wkr   = 0.0f
    var wki   = 0.0f
    var xr    = 0.0f
    var xi    = 0.0f
    var yr    = 0.0f
    var yi    = 0.0f
    var idx1  = 0
    var idx2  = 0
    m = n >> 1
    ks = 2 * nc / m
    kk = 0
    var j = 2
    while ( {
      j < m
    }) {
      k = n - j
      kk += ks
      wkr = 0.5f - c(startc + nc - kk)
      wki = c(startc + kk)
      idx1 = offa + j
      idx2 = offa + k
      xr = a(idx1) - a(idx2)
      xi = a(idx1 + 1) + a(idx2 + 1)
      yr = wkr * xr - wki * xi
      yi = wkr * xi + wki * xr
      a(idx1) -= yr
      a(idx1 + 1) -= yi
      a(idx2) += yr
      a(idx2 + 1) -= yi

      j += 2
    }
  }

  def dctsub(n: Int, a: Array[Float], offa: Int, nc: Int, c: Array[Float], startc: Int): Unit = {
    var k     = 0
    var kk    = 0
    var ks    = 0
    var m     = 0
    var wkr   = 0.0f
    var wki   = 0.0f
    var xr    = 0.0f
    var idx0  = 0
    var idx1  = 0
    var idx2  = 0
    m = n >> 1
    ks = nc / n
    kk = 0
    for (j <- 1 until m) {
      k = n - j
      kk += ks
      idx0 = startc + kk
      idx1 = offa + j
      idx2 = offa + k
      wkr = c(idx0) - c(startc + nc - kk)
      wki = c(idx0) + c(startc + nc - kk)
      xr = wki * a(idx1) - wkr * a(idx2)
      a(idx1) = wkr * a(idx1) + wki * a(idx2)
      a(idx2) = xr
    }
    a(offa + m) *= c(startc)
  }

  def scale(n: Int, m: Double, a: Array[Double], offa: Int, complex: Boolean): Unit = {
    var nthreads = ConcurrencyUtils.numThreads
    var n2 = 0
    if (complex) n2 = 2 * n
    else n2 = n
    if ((nthreads > 1) && (n2 > CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
      nthreads = 2
      val k = n2 / nthreads
      val futures = new Array[Future[_]](nthreads)
      for (i <- 0 until nthreads) {
        val firstIdx  = offa + i * k
        val lastIdx   = if (i == (nthreads - 1)) offa + n2
        else firstIdx + k
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            a(i) *= m
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      val firstIdx  = offa
      val lastIdx   = offa + n2
      for (i <- firstIdx until lastIdx) {
        a(i) *= m
      }
    }
  }

  def scale(n: Int, m: Float, a: Array[Float], offa: Int, complex: Boolean): Unit = {
    var nthreads = ConcurrencyUtils.numThreads
    var n2 = 0
    if (complex) n2 = 2 * n
    else n2 = n
    if ((nthreads > 1) && (n2 > CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
      nthreads = 2
      val k = n2 / nthreads
      val futures = new Array[Future[_]](nthreads)
      for (i <- 0 until nthreads) {
        val firstIdx  = offa + i * k
        val lastIdx   = if (i == (nthreads - 1)) offa + n2
        else firstIdx + k
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            a(i) *= m
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      val firstIdx  = offa
      val lastIdx   = offa + n2
      for (i <- firstIdx until lastIdx) {
        a(i) *= m
      }
    }
  }
}
