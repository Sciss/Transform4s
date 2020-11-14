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

package de.sciss.transform4s.fft

import de.sciss.transform4s.utils.ConcurrencyUtils

import scala.concurrent.Future
//import org.apache.commons.math3.util.FastMath._
import de.sciss.transform4s.utils.CommonUtils
import ConcurrencyUtils.executionContext

import Math.{sin, cos, ceil, log}

/**
 * Computes 1D Discrete Fourier Transform (DFT) of complex and real, double
 * precision data. The size of the data can be an arbitrary number. This is a
 * parallel implementation of split-radix and mixed-radix algorithms optimized
 * for SMP systems. <br>
 * <br>
 * This code is derived from General Purpose FFT Package written by Takuya Ooura
 * (http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html) and from JFFTPack written
 * by Baoshe Zhang (http://jfftpack.sourceforge.net/)
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
object DoubleFFT_1D {
  private final val SPLIT_RADIX = 0
  private final val MIXED_RADIX = 1
  private final val BLUESTEIN   = 2

  private final val factors     = Array(4, 2, 3, 5)
  private final val PI          = 3.14159265358979311599796346854418516
  private final val TWO_PI      = 6.28318530717958623199592693708837032

  def apply(n: Int): DoubleFFT_1D = {
    if (n < 1) throw new IllegalArgumentException("n must be greater than 0")

    var nBluestein: Int           = 0
    var ip        : Array[Int]    = null
    var w         : Array[Double] = null
    var nw        : Int           = 0
    var nc        : Int           = 0
    var wtable    : Array[Double] = null
    var wtable_r  : Array[Double] = null
    var bk1       : Array[Double] = null
    var bk2       : Array[Double] = null
    var plan      : Int           = 0

    if (!CommonUtils.isPowerOf2(n)) {
      if (CommonUtils.getReminder(n, factors) >= 211) {
        plan        = BLUESTEIN
        nBluestein  = CommonUtils.nextPow2(n * 2 - 1)
        bk1         = new Array[Double](2 * nBluestein)
        bk2         = new Array[Double](2 * nBluestein)
        ip          = new Array[Int](2 + ceil(2 + (1 << (log(nBluestein + 0.5) / log(2)).toInt / 2)).toInt)
        w           = new Array[Double](nBluestein)
        val twon    = 2 * nBluestein
        nw          = twon >> 2
        CommonUtils.makewt(nw, ip, w)
        nc          = nBluestein >> 2
        CommonUtils.makect(nc, w, nw, ip)

      } else {
        plan        = MIXED_RADIX
        wtable      = new Array[Double](4 * n + 15)
        wtable_r    = new Array[Double](2 * n + 15)
      }
    } else {
      plan          = SPLIT_RADIX
      ip            = new Array[Int](2 + ceil(2 + (1 << (log(n + 0.5) / log(2)).toInt / 2)).toInt)
      w             = new Array[Double](n)
      val twon      = 2 * n
      nw            = twon >> 2
      CommonUtils.makewt(nw, ip, w)
      nc            = n >> 2
      CommonUtils.makect(nc, w, nw, ip)
    }

    new DoubleFFT_1D(
      n           = n         ,
      nBluestein  = nBluestein,
      ip          = ip        ,
      w           = w         ,
      nw          = nw        ,
      nc          = nc        ,
      wtable      = wtable    ,
      wtable_r    = wtable_r  ,
      bk1         = bk1       ,
      bk2         = bk2       ,
      plan        = plan      ,
    )
  }
}

final class DoubleFFT_1D private (
                                   n          : Int,
                                   nBluestein : Int,
                                   ip         : Array[Int],
                                   w          : Array[Double],
                                   nw         : Int,
                                   nc         : Int,
                                   wtable     : Array[Double],
                                   wtable_r   : Array[Double],
                                   bk1        : Array[Double],
                                   bk2        : Array[Double],
                                   plan       : Int,
                        ) {
  import DoubleFFT_1D._

  // ---- init ----
  if (plan == BLUESTEIN) {
    bluesteini()
  } else if (plan == MIXED_RADIX) {
    cffti()
    rffti()
  }

  /**
   * Computes 1D forward DFT of complex data leaving the result in
   * <code>a</code>. Complex number is stored as two double values in
   * sequence: the real and imaginary part, i.e. the size of the input array
   * must be greater or equal 2*n. The physical layout of the input data has
   * to be as follows:<br>
   *
   * <pre>
   * a[2*k] = Re[k],
   * a[2*k+1] = Im[k], 0&lt;=k&lt;n
   * </pre>
   *
   * @param a data to transform
   */
  def complexForward(a: Array[Double]): Unit =
    complexForward(a, 0)

  /**
   * Computes 1D forward DFT of complex data leaving the result in
   * <code>a</code>. Complex number is stored as two double values in
   * sequence: the real and imaginary part, i.e. the size of the input array
   * must be greater or equal 2*n. The physical layout of the input data has
   * to be as follows:<br>
   *
   * <pre>
   * a[offa+2*k] = Re[k],
   * a[offa+2*k+1] = Im[k], 0&lt;=k&lt;n
   * </pre>
   *
   * @param a    data to transform
   * @param offa index of the first element in array <code>a</code>
   */
  def complexForward(a: Array[Double], offa: Int): Unit = {
    if (n == 1) return
    plan match {
      case SPLIT_RADIX =>
        CommonUtils.cftbsub(2 * n, a, offa, ip, nw, w)

      case MIXED_RADIX =>
        cfftf(a, offa, -1)

      case BLUESTEIN =>
        bluestein_complex(a, offa, -1)
    }
  }

  /**
   * Computes 1D inverse DFT of complex data leaving the result in
   * <code>a</code>. Complex number is stored as two double values in
   * sequence: the real and imaginary part, i.e. the size of the input array
   * must be greater or equal 2*n. The physical layout of the input data has
   * to be as follows:<br>
   *
   * <pre>
   * a[2*k] = Re[k],
   * a[2*k+1] = Im[k], 0&lt;=k&lt;n
   * </pre>
   *
   * @param a     data to transform
   * @param scale if true then scaling is performed
   */
  def complexInverse(a: Array[Double], scale: Boolean): Unit =
    complexInverse(a, 0, scale)

  /**
   * Computes 1D inverse DFT of complex data leaving the result in
   * <code>a</code>. Complex number is stored as two double values in
   * sequence: the real and imaginary part, i.e. the size of the input array
   * must be greater or equal 2*n. The physical layout of the input data has
   * to be as follows:<br>
   *
   * <pre>
   * a[offa+2*k] = Re[k],
   * a[offa+2*k+1] = Im[k], 0&lt;=k&lt;n
   * </pre>
   *
   * @param a     data to transform
   * @param offa  index of the first element in array <code>a</code>
   * @param scale if true then scaling is performed
   */
  def complexInverse(a: Array[Double], offa: Int, scale: Boolean): Unit = {
    if (n == 1) return
    plan match {
      case SPLIT_RADIX =>
        CommonUtils.cftfsub(2 * n, a, offa, ip, nw, w)

      case MIXED_RADIX =>
        cfftf(a, offa, +1)

      case BLUESTEIN =>
        bluestein_complex(a, offa, 1)

    }
    if (scale) CommonUtils.scale(n, 1.0 / n.toDouble, a, offa, true)
  }

  /**
   * Computes 1D forward DFT of real data leaving the result in <code>a</code>
   * . The physical layout of the output data is as follows:<br>
   *
   * if n is even then
   *
   * <pre>
   * a[2*k] = Re[k], 0&lt;=k&lt;n/2
   * a[2*k+1] = Im[k], 0&lt;k&lt;n/2
   * a[1] = Re[n/2]
   * </pre>
   *
   * if n is odd then
   *
   * <pre>
   * a[2*k] = Re[k], 0&lt;=k&lt;(n+1)/2
   * a[2*k+1] = Im[k], 0&lt;k&lt;(n-1)/2
   * a[1] = Im[(n-1)/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * forward transform, use <code>realForwardFull</code>. To get back the
   * original data, use <code>realInverse</code> on the output of this method.
   *
   * @param a data to transform
   */
  def realForward(a: Array[Double]): Unit =
    realForward(a, 0)

  /**
   * Computes 1D forward DFT of real data leaving the result in <code>a</code>
   * . The physical layout of the output data is as follows:<br>
   *
   * if n is even then
   *
   * <pre>
   * a[offa+2*k] = Re[k], 0&lt;=k&lt;n/2
   * a[offa+2*k+1] = Im[k], 0&lt;k&lt;n/2
   * a[offa+1] = Re[n/2]
   * </pre>
   *
   * if n is odd then
   *
   * <pre>
   * a[offa+2*k] = Re[k], 0&lt;=k&lt;(n+1)/2
   * a[offa+2*k+1] = Im[k], 0&lt;k&lt;(n-1)/2
   * a[offa+1] = Im[(n-1)/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * forward transform, use <code>realForwardFull</code>. To get back the
   * original data, use <code>realInverse</code> on the output of this method.
   *
   * @param a    data to transform
   * @param offa index of the first element in array <code>a</code>
   */
  def realForward(a: Array[Double], offa: Int): Unit = {
    if (n == 1) return
    plan match {
      case SPLIT_RADIX =>
        var xi = .0
        if (n > 4) {
          CommonUtils.cftfsub(n, a, offa, ip, nw, w)
          CommonUtils.rftfsub(n, a, offa, nc, w, nw)
        }
        else if (n == 4) CommonUtils.cftx020(a, offa)
        xi = a(offa) - a(offa + 1)
        a(offa) += a(offa + 1)
        a(offa + 1) = xi

      case MIXED_RADIX =>
        rfftf(a, offa)
        for (k <- n - 1 to 2 by -1) {
          val idx = offa + k
          val tmp = a(idx)
          a(idx) = a(idx - 1)
          a(idx - 1) = tmp
        }

      case BLUESTEIN =>
        bluestein_real_forward(a, offa)

    }
  }

  /**
   * Computes 1D forward DFT of real data leaving the result in <code>a</code>
   * . This method computes the full real forward transform, i.e. you will get
   * the same result as from <code>complexForward</code> called with all
   * imaginary parts equal 0. Because the result is stored in <code>a</code>,
   * the size of the input array must greater or equal 2*n, with only the
   * first n elements filled with real data. To get back the original data,
   * use <code>complexInverse</code> on the output of this method.
   *
   * @param a data to transform
   */
  def realForwardFull(a: Array[Double]): Unit =
    realForwardFull(a, 0)

  /**
   * Computes 1D forward DFT of real data leaving the result in <code>a</code>
   * . This method computes the full real forward transform, i.e. you will get
   * the same result as from <code>complexForward</code> called with all
   * imaginary part equal 0. Because the result is stored in <code>a</code>,
   * the size of the input array must greater or equal 2*n, with only the
   * first n elements filled with real data. To get back the original data,
   * use <code>complexInverse</code> on the output of this method.
   *
   * @param a    data to transform
   * @param offa index of the first element in array <code>a</code>
   */
  def realForwardFull(a: Array[Double], offa: Int): Unit = {
    val twon = 2 * n
    plan match {
      case SPLIT_RADIX =>
        realForward(a, offa)
        val nthreads = ConcurrencyUtils.numThreads
        if ((nthreads > 1) && (n / 2 > CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
          val futures = new Array[Future[_]](nthreads)
          val k = n / 2 / nthreads
          for (i <- 0 until nthreads) {
            val firstIdx = i * k
            val lastIdx = if (i == (nthreads - 1)) n / 2
            else firstIdx + k
            futures(i) = Future {
              var idx1 = 0
              var idx2 = 0
              for (k <- firstIdx until lastIdx) {
                idx1 = 2 * k
                idx2 = offa + ((twon - idx1) % twon)
                a(idx2) = a(offa + idx1)
                a(idx2 + 1) = -a(offa + idx1 + 1)
              }
            }
          }
          ConcurrencyUtils.waitForCompletion(futures)
        }
        else {
          var idx1 = 0
          var idx2 = 0
          for (k <- 0 until n / 2) {
            idx1 = 2 * k
            idx2 = offa + ((twon - idx1) % twon)
            a(idx2) = a(offa + idx1)
            a(idx2 + 1) = -a(offa + idx1 + 1)
          }
        }
        a(offa + n) = -a(offa + 1)
        a(offa + 1) = 0

      case MIXED_RADIX =>
        rfftf(a, offa)
        var m = 0
        if (n % 2 == 0) m = n / 2
        else m = (n + 1) / 2
        for (k <- 1 until m) {
          val idx1 = offa + twon - 2 * k
          val idx2 = offa + 2 * k
          a(idx1 + 1) = -a(idx2)
          a(idx1) = a(idx2 - 1)
        }
        for (k <- 1 until n) {
          val idx = offa + n - k
          val tmp = a(idx + 1)
          a(idx + 1) = a(idx)
          a(idx) = tmp
        }
        a(offa + 1) = 0

      case BLUESTEIN =>
        bluestein_real_full(a, offa, -1)

    }
  }

  /**
   * Computes 1D inverse DFT of real data leaving the result in <code>a</code>
   * . The physical layout of the input data has to be as follows:<br>
   *
   * if n is even then
   *
   * <pre>
   * a[2*k] = Re[k], 0&lt;=k&lt;n/2
   * a[2*k+1] = Im[k], 0&lt;k&lt;n/2
   * a[1] = Re[n/2]
   * </pre>
   *
   * if n is odd then
   *
   * <pre>
   * a[2*k] = Re[k], 0&lt;=k&lt;(n+1)/2
   * a[2*k+1] = Im[k], 0&lt;k&lt;(n-1)/2
   * a[1] = Im[(n-1)/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * inverse transform, use <code>realInverseFull</code>.
   *
   * @param a     data to transform
   * @param scale if true then scaling is performed
   */
  def realInverse(a: Array[Double], scale: Boolean): Unit =
    realInverse(a, 0, scale)

  /**
   * Computes 1D inverse DFT of real data leaving the result in <code>a</code>
   * . The physical layout of the input data has to be as follows:<br>
   *
   * if n is even then
   *
   * <pre>
   * a[offa+2*k] = Re[k], 0&lt;=k&lt;n/2
   * a[offa+2*k+1] = Im[k], 0&lt;k&lt;n/2
   * a[offa+1] = Re[n/2]
   * </pre>
   *
   * if n is odd then
   *
   * <pre>
   * a[offa+2*k] = Re[k], 0&lt;=k&lt;(n+1)/2
   * a[offa+2*k+1] = Im[k], 0&lt;k&lt;(n-1)/2
   * a[offa+1] = Im[(n-1)/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * inverse transform, use <code>realInverseFull</code>.
   *
   * @param a     data to transform
   * @param offa  index of the first element in array <code>a</code>
   * @param scale if true then scaling is performed
   */
  def realInverse(a: Array[Double], offa: Int, scale: Boolean): Unit = {
    if (n == 1) return
    plan match {
      case SPLIT_RADIX =>
        a(offa + 1) = 0.5 * (a(offa) - a(offa + 1))
        a(offa) -= a(offa + 1)
        if (n > 4) {
          CommonUtils.rftfsub(n, a, offa, nc, w, nw)
          CommonUtils.cftbsub(n, a, offa, ip, nw, w)
        }
        else if (n == 4) CommonUtils.cftxc020(a, offa)
        if (scale) CommonUtils.scale(n, 1.0 / (n / 2.0), a, offa, false)

      case MIXED_RADIX =>
        for (k <- 2 until n) {
          val idx = offa + k
          val tmp = a(idx - 1)
          a(idx - 1) = a(idx)
          a(idx) = tmp
        }
        rfftb(a, offa)
        if (scale) CommonUtils.scale(n, 1.0 / n, a, offa, false)

      case BLUESTEIN =>
        bluestein_real_inverse(a, offa)
        if (scale) CommonUtils.scale(n, 1.0 / n, a, offa, false)

    }
  }

  /**
   * Computes 1D inverse DFT of real data leaving the result in <code>a</code>
   * . This method computes the full real inverse transform, i.e. you will get
   * the same result as from <code>complexInverse</code> called with all
   * imaginary part equal 0. Because the result is stored in <code>a</code>,
   * the size of the input array must greater or equal 2*n, with only the
   * first n elements filled with real data.
   *
   * @param a     data to transform
   * @param scale if true then scaling is performed
   */
  def realInverseFull(a: Array[Double], scale: Boolean): Unit =
    realInverseFull(a, 0, scale)

  /**
   * Computes 1D inverse DFT of real data leaving the result in <code>a</code>
   * . This method computes the full real inverse transform, i.e. you will get
   * the same result as from <code>complexInverse</code> called with all
   * imaginary part equal 0. Because the result is stored in <code>a</code>,
   * the size of the input array must greater or equal 2*n, with only the
   * first n elements filled with real data.
   *
   * @param a     data to transform
   * @param offa  index of the first element in array <code>a</code>
   * @param scale if true then scaling is performed
   */
  def realInverseFull(a: Array[Double], offa: Int, scale: Boolean): Unit = {
    val twon = 2 * n
    plan match {
      case SPLIT_RADIX =>
        realInverse2(a, offa, scale)
        val nthreads = ConcurrencyUtils.numThreads
        if ((nthreads > 1) && (n / 2 > CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
          val futures = new Array[Future[_]](nthreads)
          val k = n / 2 / nthreads
          for (i <- 0 until nthreads) {
            val firstIdx = i * k
            val lastIdx = if (i == (nthreads - 1)) n / 2
            else firstIdx + k
            futures(i) = Future {
              var idx1 = 0
              var idx2 = 0
              for (k <- firstIdx until lastIdx) {
                idx1 = 2 * k
                idx2 = offa + ((twon - idx1) % twon)
                a(idx2) = a(offa + idx1)
                a(idx2 + 1) = -a(offa + idx1 + 1)
              }
            }
          }
          ConcurrencyUtils.waitForCompletion(futures)
        }
        else {
          var idx1 = 0
          var idx2 = 0
          for (k <- 0 until n / 2) {
            idx1 = 2 * k
            idx2 = offa + ((twon - idx1) % twon)
            a(idx2) = a(offa + idx1)
            a(idx2 + 1) = -a(offa + idx1 + 1)
          }
        }
        a(offa + n) = -a(offa + 1)
        a(offa + 1) = 0

      case MIXED_RADIX =>
        rfftf(a, offa)
        if (scale) CommonUtils.scale(n, 1.0 / n, a, offa, false)
        var m = 0
        if (n % 2 == 0) m = n / 2
        else m = (n + 1) / 2
        for (k <- 1 until m) {
          val idx1 = offa + 2 * k
          val idx2 = offa + twon - 2 * k
          a(idx1) = -a(idx1)
          a(idx2 + 1) = -a(idx1)
          a(idx2) = a(idx1 - 1)
        }
        for (k <- 1 until n) {
          val idx = offa + n - k
          val tmp = a(idx + 1)
          a(idx + 1) = a(idx)
          a(idx) = tmp
        }
        a(offa + 1) = 0

      case BLUESTEIN =>
        bluestein_real_full(a, offa, 1)
        if (scale) CommonUtils.scale(n, 1.0 / n, a, offa, true)

    }
  }

  private[fft] def realInverse2(a: Array[Double], offa: Int, scale: Boolean): Unit = {
    if (n == 1) return
    plan match {
      case SPLIT_RADIX =>
        var xi = .0
        if (n > 4) {
          CommonUtils.cftfsub(n, a, offa, ip, nw, w)
          CommonUtils.rftbsub(n, a, offa, nc, w, nw)
        }
        else if (n == 4) CommonUtils.cftbsub(n, a, offa, ip, nw, w)
        xi = a(offa) - a(offa + 1)
        a(offa) += a(offa + 1)
        a(offa + 1) = xi
        if (scale) CommonUtils.scale(n, 1.0 / n, a, offa, false)

      case MIXED_RADIX =>
        rfftf(a, offa)
        for (k <- n - 1 to 2 by -1) {
          val idx = offa + k
          val tmp = a(idx)
          a(idx) = a(idx - 1)
          a(idx - 1) = tmp
        }
        if (scale) CommonUtils.scale(n, 1.0 / n, a, offa, false)
        var m = 0
        if (n % 2 == 0) {
          m = n / 2
          for (i <- 1 until m) {
            val idx = offa + 2 * i + 1
            a(idx) = -a(idx)
          }
        }
        else {
          m = (n - 1) / 2
          for (i <- 0 until m) {
            val idx = offa + 2 * i + 1
            a(idx) = -a(idx)
          }
        }

      case BLUESTEIN =>
        bluestein_real_inverse2(a, offa)
        if (scale) CommonUtils.scale(n, 1.0 / n, a, offa, false)

    }
  }

  /* -------- initializing routines -------- */

  /*---------------------------------------------------------
     cffti: initialization of Complex FFT
     --------------------------------------------------------*/

  private[fft] def cffti(n: Int, offw: Int): Unit = {
    if (n == 1) return

    val twon  : Int     = 2 * n
    val fourn : Int     = 4 * n
    var argh  : Double  = 0.0
    var idot  : Int     = 0
    var ntry  : Int     = 0
    var i     : Int     = 0
    var j     : Int     = 0
    var argld : Double  = 0.0
    var i1    : Int     = 0
    var k1    : Int     = 0
    var l1    : Int     = 0
    var l2    : Int     = 0
    var ib    : Int     = 0
    var fi    : Double  = 0.0
    var ld    : Int     = 0
    var ii    : Int     = 0
    var nf    : Int     = 0
    var ipll  : Int     = 0
    var nll   : Int     = n
    var nq    : Int     = 0
    var nr    : Int     = 0
    var arg   : Double  = 0.0
    var ido   : Int     = 0
    var ipm   : Int     = 0

    j = 0
    // factorize_loop
    while ({
      var _continue_factorize_loop = false
      j += 1
      if (j <= 4) {
        ntry = factors(j - 1)
      } else {
        ntry += 2
      }

      while ({
        nq = nll / ntry
        nr = nll - ntry * nq
        if (nr != 0) {
          _continue_factorize_loop = true
        } else {
          nf += 1
          wtable(offw + nf + 1 + fourn) = ntry
          nll = nq
          if (ntry == 2 && nf != 1) {
            i = 2
            while (i <= nf) {
              ib = nf - i + 2
              val idx: Int = ib + fourn
              wtable(offw + idx + 1) = wtable(offw + idx)

              i += 1
            }
            wtable(offw + 2 + fourn) = 2
          }
        }

        !_continue_factorize_loop && nll != 1
      }) ()

      _continue_factorize_loop
    }) ()

    wtable(offw +     fourn) = n
    wtable(offw + 1 + fourn) = nf
    argh  = TWO_PI / n.toDouble
    i     = 1
    l1    = 1
    k1    = 1
    while (k1 <= nf) {
      ipll  = wtable(offw + k1 + 1 + fourn).toInt
      ld    = 0
      l2    = l1 * ipll
      ido   = n / l2
      idot  = ido + ido + 2
      ipm   = ipll - 1

      j = 1
      while (j <= ipm) {
        i1 = i
        wtable(offw + i - 1 + twon) = 1
        wtable(offw + i     + twon) = 0
        ld += l1
        fi = 0
        argld = ld * argh

        ii = 4
        while (ii <= idot) {
          i += 2
          fi += 1
          arg = fi * argld
          val idx: Int = i + twon
          wtable(offw + idx - 1 ) = cos(arg)
          wtable(offw + idx     ) = sin(arg)

          ii += 2
        }
        if (ipll > 5) {
          val idx1: Int = i1  + twon
          val idx2: Int = i   + twon
          wtable(offw + idx1 - 1) = wtable(offw + idx2 - 1)
          wtable(offw + idx1    ) = wtable(offw + idx2)
        }

        j += 1
      }
      l1 = l2

      k1 += 1
    }
  }

  private[fft] def cffti(): Unit = {
    if (n == 1) return

    val twon  : Int     = 2 * n
    val fourn : Int     = 4 * n
    var argh  : Double  = 0.0
    var idot  : Int     = 0
    var ntry  : Int     = 0
    var i     : Int     = 0
    var j     : Int     = 0
    var argld : Double  = 0.0
    var i1    : Int     = 0
    var k1    : Int     = 0
    var l1    : Int     = 0
    var l2    : Int     = 0
    var ib    : Int     = 0
    var fi    : Double  = 0.0
    var ld    : Int     = 0
    var ii    : Int     = 0
    var nf    : Int     = 0
    var ipll  : Int     = 0
    var nll   : Int     = n
    var nq    : Int     = 0
    var nr    : Int     = 0
    var arg   : Double  = 0.0
    var ido   : Int     = 0
    var ipm   : Int     = 0

    j = 0
    // factorize_loop
    while ({
      var _continue_factorize_loop = false
      j += 1
      if (j <= 4) {
        ntry = factors(j - 1)
      } else {
        ntry += 2
      }

      while ({
        nq = nll / ntry
        nr = nll - ntry * nq
        if (nr != 0) {
          _continue_factorize_loop = true
        } else {
          nf += 1
          wtable(nf + 1 + fourn) = ntry
          nll = nq
          if (ntry == 2 && nf != 1) {
            i = 2
            while (i <= nf) {
              ib = nf - i + 2
              val idx: Int = ib + fourn
              wtable(idx + 1) = wtable(idx)

              i += 1
            }
            wtable(2 + fourn) = 2
          }
        }

        !_continue_factorize_loop && nll != 1
      }) ()

      _continue_factorize_loop
    }) ()

    wtable(fourn)     = n
    wtable(1 + fourn) = nf
    argh = TWO_PI / n.toDouble
    i   = 1
    l1  = 1
    k1  = 1
    while (k1 <= nf) {
      ipll  = wtable(k1 + 1 + fourn).toInt
      ld    = 0
      l2    = l1 * ipll
      ido   = n / l2
      idot  = ido + ido + 2
      ipm   = ipll - 1

      j = 1
      while (j <= ipm) {
        i1 = i
        wtable(i - 1  + twon) = 1
        wtable(i      + twon) = 0
        ld += l1
        fi = 0
        argld = ld * argh

        ii = 4
        while (ii <= idot) {
          i  += 2
          fi += 1
          arg = fi * argld
          val idx: Int = i + twon
          wtable(idx - 1) = cos(arg)
          wtable(idx    ) = sin(arg)

          ii += 2
        }
        if (ipll > 5) {
          val idx1: Int = i1  + twon
          val idx2: Int = i   + twon
          wtable(idx1 - 1 ) = wtable(idx2 - 1)
          wtable(idx1     ) = wtable(idx2)
        }

        j += 1
      }
      l1 = l2

      k1 += 1
    }
  }

  private[fft] def rffti(): Unit = {
    if (n == 1) return

    val twon  : Int     = 2 * n
    var argh  : Double  = 0.0
    var ntry  : Int     = 0
    var i     : Int     = 0
    var j     : Int     = 0
    var argld : Double  = 0.0
    var k1    : Int     = 0
    var l1    : Int     = 0
    var l2    : Int     = 0
    var ib    : Int     = 0
    var fi    : Double  = 0.0
    var ld    : Int     = 0
    var ii    : Int     = 0
    var nf    : Int     = 0
    var ipll  : Int     = 0
    var nll   : Int     = n
    var is    : Int     = 0
    var nq    : Int     = 0
    var nr    : Int     = 0
    var arg   : Double  = 0.0
    var ido   : Int     = 0
    var ipm   : Int     = 0
    var nfm1  : Int     = 0

    j = 0
    // factorize_loop
    while ({
      var _continue_factorize_loop = false
      j += 1
      if (j <= 4) {
        ntry = factors(j - 1)
      } else {
        ntry += 2
      }

      while ({
        nq = nll / ntry
        nr = nll - ntry * nq
        if (nr != 0) {
          _continue_factorize_loop = true
        } else {
          nf += 1
          wtable_r(nf + 1 + twon) = ntry
          nll = nq
          if (ntry == 2 && nf != 1) {
            i = 2
            while (i <= nf) {
              ib = nf - i + 2
              val idx: Int = ib + twon
              wtable_r(idx + 1) = wtable_r(idx)

              i += 1
            }
            wtable_r(2 + twon) = 2
          }
        }

        !_continue_factorize_loop && nll != 1
      }) ()

      _continue_factorize_loop
    }) ()

    wtable_r(twon) = n
    wtable_r(1 + twon) = nf
    argh = TWO_PI / n.toDouble
    is = 0
    nfm1 = nf - 1
    l1 = 1
    if (nfm1 == 0) {
      return
    }
    k1 = 1
    while ( {
      k1 <= nfm1
    }) {
      ipll = wtable_r(k1 + 1 + twon).toInt
      ld = 0
      l2 = l1 * ipll
      ido = n / l2
      ipm = ipll - 1
      j = 1
      while ( {
        j <= ipm
      }) {
        ld += l1
        i = is
        argld = ld.toDouble * argh
        fi = 0
        ii = 3
        while ( {
          ii <= ido
        }) {
          i += 2
          fi += 1
          arg = fi * argld
          val idx: Int = i + n
          wtable_r(idx - 2) = cos(arg)
          wtable_r(idx - 1) = sin(arg)

          ii += 2
        }
        is += ido

        j += 1
      }
      l1 = l2

      k1 += 1
    }
  }

  private def bluesteini(): Unit = {
    var k: Int = 0
    var arg: Double = .0
    val pi_n: Double = PI / n
    bk1(0) = 1
    bk1(1) = 0
    for (i <- 1 until n) {
      k += 2 * i - 1
      if (k >= 2 * n) {
        k -= 2 * n
      }
      arg = pi_n * k
      bk1(2 * i     ) = cos(arg)
      bk1(2 * i + 1 ) = sin(arg)
    }
    val scale: Double = 1.0 / nBluestein
    bk2(0) = bk1(0) * scale
    bk2(1) = bk1(1) * scale
    var i: Int = 2
    while ( {
      i < 2 * n
    }) {
      bk2(i) = bk1(i) * scale
      bk2(i + 1) = bk1(i + 1) * scale
      bk2(2 * nBluestein - i) = bk2(i)
      bk2(2 * nBluestein - i + 1) = bk2(i + 1)

      i += 2
    }
    CommonUtils.cftbsub(2 * nBluestein, bk2, 0, ip, nw, w)
  }

  private def bluestein_complex(a: Array[Double], offa: Int, isign: Int): Unit = {
    val ak: Array[Double] = new Array[Double](2 * nBluestein)
    val threads: Int = ConcurrencyUtils.numThreads
    if ((threads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
      var nthreads: Int = 2
      if ((threads >= 4) && (n >= CommonUtils.threadsBeginN_1D_FFT_4Threads)) {
        nthreads = 4
      }
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var k: Int = n / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          n
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          if (isign > 0) {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val idx3: Int = offa + idx1
              val idx4: Int = offa + idx2
              ak(idx1) = a(idx3) * bk1(idx1) - a(idx4) * bk1(idx2)
              ak(idx2) = a(idx3) * bk1(idx2) + a(idx4) * bk1(idx1)
            }
          }
          else {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val idx3: Int = offa + idx1
              val idx4: Int = offa + idx2
              ak(idx1) = a(idx3) * bk1(idx1) + a(idx4) * bk1(idx2)
              ak(idx2) = -a(idx3) * bk1(idx2) + a(idx4) * bk1(idx1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      k = nBluestein / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          nBluestein
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          if (isign > 0) {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
              ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
              ak(idx2) = im
            }
          }
          else {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val im: Double = ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
              ak(idx1) = ak(idx1) * bk2(idx1) - ak(idx2) * bk2(idx2)
              ak(idx2) = im
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
      k = n / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          n
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          if (isign > 0) {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val idx3: Int = offa + idx1
              val idx4: Int = offa + idx2
              a(idx3) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
              a(idx4) = bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
            }
          }
          else {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val idx3: Int = offa + idx1
              val idx4: Int = offa + idx2
              a(idx3) = bk1(idx1) * ak(idx1) + bk1(idx2) * ak(idx2)
              a(idx4) = -bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      if (isign > 0) {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val idx3: Int = offa + idx1
          val idx4: Int = offa + idx2
          ak(idx1) = a(idx3) * bk1(idx1) - a(idx4) * bk1(idx2)
          ak(idx2) = a(idx3) * bk1(idx2) + a(idx4) * bk1(idx1)
        }
      }
      else {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val idx3: Int = offa + idx1
          val idx4: Int = offa + idx2
          ak(idx1) = a(idx3) * bk1(idx1) + a(idx4) * bk1(idx2)
          ak(idx2) = -a(idx3) * bk1(idx2) + a(idx4) * bk1(idx1)
        }
      }
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      if (isign > 0) {
        for (i <- 0 until nBluestein) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
          ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
          ak(idx2) = im
        }
      }
      else {
        for (i <- 0 until nBluestein) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val im: Double = ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
          ak(idx1) = ak(idx1) * bk2(idx1) - ak(idx2) * bk2(idx2)
          ak(idx2) = im
        }
      }
      CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
      if (isign > 0) {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val idx3: Int = offa + idx1
          val idx4: Int = offa + idx2
          a(idx3) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
          a(idx4) = bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
        }
      }
      else {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val idx3: Int = offa + idx1
          val idx4: Int = offa + idx2
          a(idx3) = bk1(idx1) * ak(idx1) + bk1(idx2) * ak(idx2)
          a(idx4) = -bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
        }
      }
    }
  }

  private def bluestein_real_full(a: Array[Double], offa: Int, isign: Int): Unit = {
    val ak: Array[Double] = new Array[Double](2 * nBluestein)
    val threads: Int = ConcurrencyUtils.numThreads
    if ((threads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
      var nthreads: Int = 2
      if ((threads >= 4) && (n >= CommonUtils.threadsBeginN_1D_FFT_4Threads)) {
        nthreads = 4
      }
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var k: Int = n / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          n
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          if (isign > 0) {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val idx3: Int = offa + i
              ak(idx1) = a(idx3) * bk1(idx1)
              ak(idx2) = a(idx3) * bk1(idx2)
            }
          }
          else {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val idx3: Int = offa + i
              ak(idx1) = a(idx3) * bk1(idx1)
              ak(idx2) = -a(idx3) * bk1(idx2)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      k = nBluestein / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          nBluestein
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          if (isign > 0) {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
              ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
              ak(idx2) = im
            }
          }
          else {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              val im: Double = ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
              ak(idx1) = ak(idx1) * bk2(idx1) - ak(idx2) * bk2(idx2)
              ak(idx2) = im
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
      k = n / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          n
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          if (isign > 0) {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              a(offa + idx1) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
              a(offa + idx2) = bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
            }
          }
          else {
            for (i <- firstIdx until lastIdx) {
              val idx1: Int = 2 * i
              val idx2: Int = idx1 + 1
              a(offa + idx1) = bk1(idx1) * ak(idx1) + bk1(idx2) * ak(idx2)
              a(offa + idx2) = -bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      if (isign > 0) {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val idx3: Int = offa + i
          ak(idx1) = a(idx3) * bk1(idx1)
          ak(idx2) = a(idx3) * bk1(idx2)
        }
      }
      else {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val idx3: Int = offa + i
          ak(idx1) = a(idx3) * bk1(idx1)
          ak(idx2) = -a(idx3) * bk1(idx2)
        }
      }
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      if (isign > 0) {
        for (i <- 0 until nBluestein) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
          ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
          ak(idx2) = im
        }
      }
      else {
        for (i <- 0 until nBluestein) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          val im: Double = ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
          ak(idx1) = ak(idx1) * bk2(idx1) - ak(idx2) * bk2(idx2)
          ak(idx2) = im
        }
      }
      CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
      if (isign > 0) {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          a(offa + idx1) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
          a(offa + idx2) = bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
        }
      }
      else {
        for (i <- 0 until n) {
          val idx1: Int = 2 * i
          val idx2: Int = idx1 + 1
          a(offa + idx1) = bk1(idx1) * ak(idx1) + bk1(idx2) * ak(idx2)
          a(offa + idx2) = -bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
        }
      }
    }
  }

  private def bluestein_real_forward(a: Array[Double], offa: Int): Unit = {
    val ak: Array[Double] = new Array[Double](2 * nBluestein)
    val threads: Int = ConcurrencyUtils.numThreads
    if ((threads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
      var nthreads: Int = 2
      if ((threads >= 4) && (n >= CommonUtils.threadsBeginN_1D_FFT_4Threads)) {
        nthreads = 4
      }
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var k: Int = n / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          n
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            val idx1: Int = 2 * i
            val idx2: Int = idx1 + 1
            val idx3: Int = offa + i
            ak(idx1) = a(idx3) * bk1(idx1)
            ak(idx2) = -a(idx3) * bk1(idx2)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      k = nBluestein / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          nBluestein
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            val idx1: Int = 2 * i
            val idx2: Int = idx1 + 1
            val im: Double = ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
            ak(idx1) = ak(idx1) * bk2(idx1) - ak(idx2) * bk2(idx2)
            ak(idx2) = im
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (i <- 0 until n) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val idx3: Int = offa + i
        ak(idx1) = a(idx3) * bk1(idx1)
        ak(idx2) = -a(idx3) * bk1(idx2)
      }
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      for (i <- 0 until nBluestein) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val im: Double = ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
        ak(idx1) = ak(idx1) * bk2(idx1) - ak(idx2) * bk2(idx2)
        ak(idx2) = im
      }
    }
    CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
    if (n % 2 == 0) {
      a(offa) = bk1(0) * ak(0) + bk1(1) * ak(1)
      a(offa + 1) = bk1(n) * ak(n) + bk1(n + 1) * ak(n + 1)
      for (i <- 1 until n / 2) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        a(offa + idx1) = bk1(idx1) * ak(idx1) + bk1(idx2) * ak(idx2)
        a(offa + idx2) = -bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
      }
    }
    else {
      a(offa) = bk1(0) * ak(0) + bk1(1) * ak(1)
      a(offa + 1) = -bk1(n) * ak(n - 1) + bk1(n - 1) * ak(n)
      for (i <- 1 until (n - 1) / 2) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        a(offa + idx1) = bk1(idx1) * ak(idx1) + bk1(idx2) * ak(idx2)
        a(offa + idx2) = -bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
      }
      a(offa + n - 1) = bk1(n - 1) * ak(n - 1) + bk1(n) * ak(n)
    }
  }

  private def bluestein_real_inverse(a: Array[Double], offa: Int): Unit = {
    val ak: Array[Double] = new Array[Double](2 * nBluestein)
    if (n % 2 == 0) {
      ak(0) = a(offa) * bk1(0)
      ak(1) = a(offa) * bk1(1)
      for (i <- 1 until n / 2) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val idx3: Int = offa + idx1
        val idx4: Int = offa + idx2
        ak(idx1) = a(idx3) * bk1(idx1) - a(idx4) * bk1(idx2)
        ak(idx2) = a(idx3) * bk1(idx2) + a(idx4) * bk1(idx1)
      }
      ak(n) = a(offa + 1) * bk1(n)
      ak(n + 1) = a(offa + 1) * bk1(n + 1)
      for (i <- n / 2 + 1 until n) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val idx3: Int = offa + 2 * n - idx1
        val idx4: Int = idx3 + 1
        ak(idx1) = a(idx3) * bk1(idx1) + a(idx4) * bk1(idx2)
        ak(idx2) = a(idx3) * bk1(idx2) - a(idx4) * bk1(idx1)
      }
    }
    else {
      ak(0) = a(offa) * bk1(0)
      ak(1) = a(offa) * bk1(1)
      for (i <- 1 until (n - 1) / 2) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val idx3: Int = offa + idx1
        val idx4: Int = offa + idx2
        ak(idx1) = a(idx3) * bk1(idx1) - a(idx4) * bk1(idx2)
        ak(idx2) = a(idx3) * bk1(idx2) + a(idx4) * bk1(idx1)
      }
      ak(n - 1) = a(offa + n - 1) * bk1(n - 1) - a(offa + 1) * bk1(n)
      ak(n) = a(offa + n - 1) * bk1(n) + a(offa + 1) * bk1(n - 1)
      ak(n + 1) = a(offa + n - 1) * bk1(n + 1) + a(offa + 1) * bk1(n + 2)
      ak(n + 2) = a(offa + n - 1) * bk1(n + 2) - a(offa + 1) * bk1(n + 1)
      for (i <- (n - 1) / 2 + 2 until n) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val idx3: Int = offa + 2 * n - idx1
        val idx4: Int = idx3 + 1
        ak(idx1) = a(idx3) * bk1(idx1) + a(idx4) * bk1(idx2)
        ak(idx2) = a(idx3) * bk1(idx2) - a(idx4) * bk1(idx1)
      }
    }
    CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
    val threads: Int = ConcurrencyUtils.numThreads
    if ((threads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
      var nthreads: Int = 2
      if ((threads >= 4) && (n >= CommonUtils.threadsBeginN_1D_FFT_4Threads)) {
        nthreads = 4
      }
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var k: Int = nBluestein / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          nBluestein
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            val idx1: Int = 2 * i
            val idx2: Int = idx1 + 1
            val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
            ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
            ak(idx2) = im
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
      k = n / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          n
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            val idx1: Int = 2 * i
            val idx2: Int = idx1 + 1
            a(offa + i) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (i <- 0 until nBluestein) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
        ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
        ak(idx2) = im
      }
      CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
      for (i <- 0 until n) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        a(offa + i) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
      }
    }
  }

  private def bluestein_real_inverse2(a: Array[Double], offa: Int): Unit = {
    val ak: Array[Double] = new Array[Double](2 * nBluestein)
    val threads: Int = ConcurrencyUtils.numThreads
    if ((threads > 1) && (n >= CommonUtils.threadsBeginN_1D_FFT_2Threads)) {
      var nthreads: Int = 2
      if ((threads >= 4) && (n >= CommonUtils.threadsBeginN_1D_FFT_4Threads)) {
        nthreads = 4
      }
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var k: Int = n / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          n
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            val idx1: Int = 2 * i
            val idx2: Int = idx1 + 1
            val idx3: Int = offa + i
            ak(idx1) = a(idx3) * bk1(idx1)
            ak(idx2) = a(idx3) * bk1(idx2)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      k = nBluestein / nthreads
      for (i <- 0 until nthreads) {
        val firstIdx: Int = i * k
        val lastIdx: Int = if (i == (nthreads - 1)) {
          nBluestein
        }
        else {
          firstIdx + k
        }
        futures(i) = Future {
          for (i <- firstIdx until lastIdx) {
            val idx1: Int = 2 * i
            val idx2: Int = idx1 + 1
            val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
            ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
            ak(idx2) = im
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (i <- 0 until n) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val idx3: Int = offa + i
        ak(idx1) = a(idx3) * bk1(idx1)
        ak(idx2) = a(idx3) * bk1(idx2)
      }
      CommonUtils.cftbsub(2 * nBluestein, ak, 0, ip, nw, w)
      for (i <- 0 until nBluestein) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        val im: Double = -ak(idx1) * bk2(idx2) + ak(idx2) * bk2(idx1)
        ak(idx1) = ak(idx1) * bk2(idx1) + ak(idx2) * bk2(idx2)
        ak(idx2) = im
      }
    }
    CommonUtils.cftfsub(2 * nBluestein, ak, 0, ip, nw, w)
    if (n % 2 == 0) {
      a(offa) = bk1(0) * ak(0) - bk1(1) * ak(1)
      a(offa + 1) = bk1(n) * ak(n) - bk1(n + 1) * ak(n + 1)
      for (i <- 1 until n / 2) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        a(offa + idx1) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
        a(offa + idx2) = bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
      }
    }
    else {
      a(offa) = bk1(0) * ak(0) - bk1(1) * ak(1)
      a(offa + 1) = bk1(n) * ak(n - 1) + bk1(n - 1) * ak(n)
      for (i <- 1 until (n - 1) / 2) {
        val idx1: Int = 2 * i
        val idx2: Int = idx1 + 1
        a(offa + idx1) = bk1(idx1) * ak(idx1) - bk1(idx2) * ak(idx2)
        a(offa + idx2) = bk1(idx2) * ak(idx1) + bk1(idx1) * ak(idx2)
      }
      a(offa + n - 1) = bk1(n - 1) * ak(n - 1) - bk1(n) * ak(n)
    }
  }

  /*---------------------------------------------------------
     rfftf1: further processing of Real forward FFT
     --------------------------------------------------------*/
  private[fft] def rfftf(a: Array[Double], offa: Int): Unit = {
    if (n == 1) return

    val ch    = new Array[Double](n)
    val twon  = 2 * n
    val nf    = wtable_r(1 + twon).toInt

    var l1    = 0
    var l2    = n
    var na    = 1
    var kh    = 0
    var ipll  = 0
    var iw    = twon - 1
    var ido   = 0
    var idl1  = 0

    for (k1 <- 1 to nf) {
      kh    = nf - k1
      ipll  = wtable_r(kh + 2 + twon).toInt
      l1    = l2 / ipll
      ido   = n / l2
      idl1  = ido * l1
      iw   -= (ipll - 1) * ido
      na    = 1 - na
      ipll match {
        case 2 =>
          if (na == 0)  radf2(ido, l1, a, offa, ch, 0, iw)
          else          radf2(ido, l1, ch, 0, a, offa, iw)

        case 3 =>
          if (na == 0)  radf3(ido, l1, a, offa, ch, 0, iw)
          else          radf3(ido, l1, ch, 0, a, offa, iw)

        case 4 =>
          if (na == 0)  radf4(ido, l1, a, offa, ch, 0, iw)
          else          radf4(ido, l1, ch, 0, a, offa, iw)

        case 5 =>
          if (na == 0)  radf5(ido, l1, a, offa, ch, 0, iw)
          else          radf5(ido, l1, ch, 0, a, offa, iw)

        case _ =>
          if (ido == 1) na = 1 - na
          if (na == 0) {
            radfg(ido, ipll, l1, idl1, a, offa, ch, 0, iw)
            na = 1
          }
          else {
            radfg(ido, ipll, l1, idl1, ch, 0, a, offa, iw)
            na = 0
          }

      }
      l2 = l1
    }
    if (na == 1) return
    System.arraycopy(ch, 0, a, offa, n)
  }

  /*---------------------------------------------------------
   rfftb1: further processing of Real backward FFT
   --------------------------------------------------------*/
  private[fft] def rfftb(a: Array[Double], offa: Int): Unit = {
    if (n == 1) return
    var l1 = 0
    var l2 = 0
    var na = 0
    var nf = 0
    var ipll = 0
    var iw = 0
    var ido = 0
    var idl1 = 0
    val ch = new Array[Double](n)
    val twon = 2 * n
    nf = wtable_r(1 + twon).toInt
    na = 0
    l1 = 1
    iw = n
    for (k1 <- 1 to nf) {
      ipll = wtable_r(k1 + 1 + twon).toInt
      l2 = ipll * l1
      ido = n / l2
      idl1 = ido * l1
      ipll match {
        case 2 =>
          if (na == 0) radb2(ido, l1, a, offa, ch, 0, iw)
          else radb2(ido, l1, ch, 0, a, offa, iw)
          na = 1 - na

        case 3 =>
          if (na == 0) radb3(ido, l1, a, offa, ch, 0, iw)
          else radb3(ido, l1, ch, 0, a, offa, iw)
          na = 1 - na

        case 4 =>
          if (na == 0) radb4(ido, l1, a, offa, ch, 0, iw)
          else radb4(ido, l1, ch, 0, a, offa, iw)
          na = 1 - na

        case 5 =>
          if (na == 0) radb5(ido, l1, a, offa, ch, 0, iw)
          else radb5(ido, l1, ch, 0, a, offa, iw)
          na = 1 - na

        case _ =>
          if (na == 0) radbg(ido, ipll, l1, idl1, a, offa, ch, 0, iw)
          else radbg(ido, ipll, l1, idl1, ch, 0, a, offa, iw)
          if (ido == 1) na = 1 - na

      }
      l1 = l2
      iw += (ipll - 1) * ido
    }
    if (na == 0) return
    System.arraycopy(ch, 0, a, offa, n)
  }

  /*-------------------------------------------------
   radf2: Real FFT's forward processing of factor 2
   -------------------------------------------------*/
  private[fft] def radf2(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    var i = 0
    var ic = 0
    var idx0 = 0
    var idx1 = 0
    var idx2 = 0
    var idx3 = 0
    var idx4 = 0
    var t1i = .0
    var t1r = .0
    var w1r = .0
    var w1i = .0
    var iw1 = 0
    iw1 = offset
    idx0 = l1 * ido
    idx1 = 2 * ido
    for (k <- 0 until l1) {
      val oidx1 = out_off + k * idx1
      val oidx2 = oidx1 + idx1 - 1
      val iidx1 = in_off + k * ido
      val iidx2 = iidx1 + idx0
      val i1r = in(iidx1)
      val i2r = in(iidx2)
      out(oidx1) = i1r + i2r
      out(oidx2) = i1r - i2r
    }
    if (ido < 2) return
    if (ido != 2) {
      for (k <- 0 until l1) {
        idx1 = k * ido
        idx2 = 2 * idx1
        idx3 = idx2 + ido
        idx4 = idx1 + idx0
        i = 2
        while ( {
          i < ido
        }) {
          ic = ido - i
          val widx1 = i - 1 + iw1
          val oidx1 = out_off + i + idx2
          val oidx2 = out_off + ic + idx3
          val iidx1 = in_off + i + idx1
          val iidx2 = in_off + i + idx4
          val a1i = in(iidx1 - 1)
          val a1r = in(iidx1)
          val a2i = in(iidx2 - 1)
          val a2r = in(iidx2)
          w1r = wtable_r(widx1 - 1)
          w1i = wtable_r(widx1)
          t1r = w1r * a2i + w1i * a2r
          t1i = w1r * a2r - w1i * a2i
          out(oidx1) = a1r + t1i
          out(oidx1 - 1) = a1i + t1r
          out(oidx2) = t1i - a1r
          out(oidx2 - 1) = a1i - t1r

          i += 2
        }
      }
      if (ido % 2 == 1) return
    }
    idx2 = 2 * idx1
    for (k <- 0 until l1) {
      idx1 = k * ido
      val oidx1 = out_off + idx2 + ido
      val iidx1 = in_off + ido - 1 + idx1
      out(oidx1) = -in(iidx1 + idx0)
      out(oidx1 - 1) = in(iidx1)
    }
  }

  /*-------------------------------------------------
   radb2: Real FFT's backward processing of factor 2
   -------------------------------------------------*/
  private[fft] def radb2(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    var i = 0
    var ic = 0
    var t1i = .0
    var t1r = .0
    var w1r = .0
    var w1i = .0
    val iw1 = offset
    val idx0 = l1 * ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 2 * idx1
      val idx3 = idx2 + ido
      val oidx1 = out_off + idx1
      val iidx1 = in_off + idx2
      val iidx2 = in_off + ido - 1 + idx3
      val i1r = in(iidx1)
      val i2r = in(iidx2)
      out(oidx1) = i1r + i2r
      out(oidx1 + idx0) = i1r - i2r
    }
    if (ido < 2) return
    if (ido != 2) {
      for (k <- 0 until l1) {
        val idx1 = k * ido
        val idx2 = 2 * idx1
        val idx3 = idx2 + ido
        val idx4 = idx1 + idx0
        i = 2
        while ( {
          i < ido
        }) {
          ic = ido - i
          val idx5 = i - 1 + iw1
          val idx6 = out_off + i
          val idx7 = in_off + i
          val idx8 = in_off + ic
          w1r = wtable_r(idx5 - 1)
          w1i = wtable_r(idx5)
          val iidx1 = idx7 + idx2
          val iidx2 = idx8 + idx3
          val oidx1 = idx6 + idx1
          val oidx2 = idx6 + idx4
          t1r = in(iidx1 - 1) - in(iidx2 - 1)
          t1i = in(iidx1) + in(iidx2)
          val i1i = in(iidx1)
          val i1r = in(iidx1 - 1)
          val i2i = in(iidx2)
          val i2r = in(iidx2 - 1)
          out(oidx1 - 1) = i1r + i2r
          out(oidx1) = i1i - i2i
          out(oidx2 - 1) = w1r * t1r - w1i * t1i
          out(oidx2) = w1r * t1i + w1i * t1r

          i += 2
        }
      }
      if (ido % 2 == 1) return
    }
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 2 * idx1
      val oidx1 = out_off + ido - 1 + idx1
      val iidx1 = in_off + idx2 + ido
      out(oidx1) = 2 * in(iidx1 - 1)
      out(oidx1 + idx0) = -2 * in(iidx1)
    }
  }

  /*-------------------------------------------------
   radf3: Real FFT's forward processing of factor 3
   -------------------------------------------------*/
  private[fft] def radf3(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    val taur = -0.5
    val taui = 0.866025403784438707610604524234076962
    var i = 0
    var ic = 0
    var ci2 = .0
    var di2 = .0
    var di3 = .0
    var cr2 = .0
    var dr2 = .0
    var dr3 = .0
    var ti2 = .0
    var ti3 = .0
    var tr2 = .0
    var tr3 = .0
    var w1r = .0
    var w2r = .0
    var w1i = .0
    var w2i = .0
    var iw1 = 0
    var iw2 = 0
    iw1 = offset
    iw2 = iw1 + ido
    val idx0 = l1 * ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx3 = 2 * idx0
      val idx4 = (3 * k + 1) * ido
      val iidx1 = in_off + idx1
      val iidx2 = iidx1 + idx0
      val iidx3 = iidx1 + idx3
      val i1r = in(iidx1)
      val i2r = in(iidx2)
      val i3r = in(iidx3)
      cr2 = i2r + i3r
      out(out_off + 3 * idx1) = i1r + cr2
      out(out_off + idx4 + ido) = taui * (i3r - i2r)
      out(out_off + ido - 1 + idx4) = i1r + taur * cr2
    }
    if (ido == 1) return
    for (k <- 0 until l1) {
      val idx3 = k * ido
      val idx4 = 3 * idx3
      val idx5 = idx3 + idx0
      val idx6 = idx5 + idx0
      val idx7 = idx4 + ido
      val idx8 = idx7 + ido
      i = 2
      while ( {
        i < ido
      }) {
        ic = ido - i
        val widx1 = i - 1 + iw1
        val widx2 = i - 1 + iw2
        w1r = wtable_r(widx1 - 1)
        w1i = wtable_r(widx1)
        w2r = wtable_r(widx2 - 1)
        w2i = wtable_r(widx2)
        val idx9 = in_off + i
        val idx10 = out_off + i
        val idx11 = out_off + ic
        val iidx1 = idx9 + idx3
        val iidx2 = idx9 + idx5
        val iidx3 = idx9 + idx6
        val i1i = in(iidx1 - 1)
        val i1r = in(iidx1)
        val i2i = in(iidx2 - 1)
        val i2r = in(iidx2)
        val i3i = in(iidx3 - 1)
        val i3r = in(iidx3)
        dr2 = w1r * i2i + w1i * i2r
        di2 = w1r * i2r - w1i * i2i
        dr3 = w2r * i3i + w2i * i3r
        di3 = w2r * i3r - w2i * i3i
        cr2 = dr2 + dr3
        ci2 = di2 + di3
        tr2 = i1i + taur * cr2
        ti2 = i1r + taur * ci2
        tr3 = taui * (di2 - di3)
        ti3 = taui * (dr3 - dr2)
        val oidx1 = idx10 + idx4
        val oidx2 = idx11 + idx7
        val oidx3 = idx10 + idx8
        out(oidx1 - 1) = i1i + cr2
        out(oidx1) = i1r + ci2
        out(oidx2 - 1) = tr2 - tr3
        out(oidx2) = ti3 - ti2
        out(oidx3 - 1) = tr2 + tr3
        out(oidx3) = ti2 + ti3

        i += 2
      }
    }
  }

  /*-------------------------------------------------
     radb3: Real FFT's backward processing of factor 3
     -------------------------------------------------*/
  private[fft] def radb3(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    val taur = -0.5
    val taui = 0.866025403784438707610604524234076962
    var i = 0
    var ic = 0
    var ci2 = .0
    var ci3 = .0
    var di2 = .0
    var di3 = .0
    var cr2 = .0
    var cr3 = .0
    var dr2 = .0
    var dr3 = .0
    var ti2 = .0
    var tr2 = .0
    var w1r = .0
    var w2r = .0
    var w1i = .0
    var w2i = .0
    var iw1 = 0
    var iw2 = 0
    iw1 = offset
    iw2 = iw1 + ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val iidx1 = in_off + 3 * idx1
      val iidx2 = iidx1 + 2 * ido
      val i1i = in(iidx1)
      tr2 = 2 * in(iidx2 - 1)
      cr2 = i1i + taur * tr2
      ci3 = 2 * taui * in(iidx2)
      out(out_off + idx1) = i1i + tr2
      out(out_off + (k + l1) * ido) = cr2 - ci3
      out(out_off + (k + 2 * l1) * ido) = cr2 + ci3
    }
    if (ido == 1) return
    val idx0 = l1 * ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 3 * idx1
      val idx3 = idx2 + ido
      val idx4 = idx3 + ido
      val idx5 = idx1 + idx0
      val idx6 = idx5 + idx0
      i = 2
      while ( {
        i < ido
      }) {
        ic = ido - i
        val idx7 = in_off + i
        val idx8 = in_off + ic
        val idx9 = out_off + i
        val iidx1 = idx7 + idx2
        val iidx2 = idx7 + idx4
        val iidx3 = idx8 + idx3
        val i1i = in(iidx1 - 1)
        val i1r = in(iidx1)
        val i2i = in(iidx2 - 1)
        val i2r = in(iidx2)
        val i3i = in(iidx3 - 1)
        val i3r = in(iidx3)
        tr2 = i2i + i3i
        cr2 = i1i + taur * tr2
        ti2 = i2r - i3r
        ci2 = i1r + taur * ti2
        cr3 = taui * (i2i - i3i)
        ci3 = taui * (i2r + i3r)
        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3
        val widx1 = i - 1 + iw1
        val widx2 = i - 1 + iw2
        w1r = wtable_r(widx1 - 1)
        w1i = wtable_r(widx1)
        w2r = wtable_r(widx2 - 1)
        w2i = wtable_r(widx2)
        val oidx1 = idx9 + idx1
        val oidx2 = idx9 + idx5
        val oidx3 = idx9 + idx6
        out(oidx1 - 1) = i1i + tr2
        out(oidx1) = i1r + ti2
        out(oidx2 - 1) = w1r * dr2 - w1i * di2
        out(oidx2) = w1r * di2 + w1i * dr2
        out(oidx3 - 1) = w2r * dr3 - w2i * di3
        out(oidx3) = w2r * di3 + w2i * dr3

        i += 2
      }
    }
  }

  /*-------------------------------------------------
   radf4: Real FFT's forward processing of factor 4
   -------------------------------------------------*/
  private[fft] def radf4(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    val hsqt2 = 0.707106781186547572737310929369414225
    var i = 0
    var ic = 0
    var ci2 = .0
    var ci3 = .0
    var ci4 = .0
    var cr2 = .0
    var cr3 = .0
    var cr4 = .0
    var ti1 = .0
    var ti2 = .0
    var ti3 = .0
    var ti4 = .0
    var tr1 = .0
    var tr2 = .0
    var tr3 = .0
    var tr4 = .0
    var w1r = .0
    var w1i = .0
    var w2r = .0
    var w2i = .0
    var w3r = .0
    var w3i = .0
    var iw1 = 0
    var iw2 = 0
    var iw3 = 0
    iw1 = offset
    iw2 = offset + ido
    iw3 = iw2 + ido
    val idx0 = l1 * ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 4 * idx1
      val idx3 = idx1 + idx0
      val idx4 = idx3 + idx0
      val idx5 = idx4 + idx0
      val idx6 = idx2 + ido
      val i1r = in(in_off + idx1)
      val i2r = in(in_off + idx3)
      val i3r = in(in_off + idx4)
      val i4r = in(in_off + idx5)
      tr1 = i2r + i4r
      tr2 = i1r + i3r
      val oidx1 = out_off + idx2
      val oidx2 = out_off + idx6 + ido
      out(oidx1) = tr1 + tr2
      out(oidx2 - 1 + ido + ido) = tr2 - tr1
      out(oidx2 - 1) = i1r - i3r
      out(oidx2) = i4r - i2r
    }
    if (ido < 2) return
    if (ido != 2) {
      for (k <- 0 until l1) {
        val idx1 = k * ido
        val idx2 = idx1 + idx0
        val idx3 = idx2 + idx0
        val idx4 = idx3 + idx0
        val idx5 = 4 * idx1
        val idx6 = idx5 + ido
        val idx7 = idx6 + ido
        val idx8 = idx7 + ido
        i = 2
        while ( {
          i < ido
        }) {
          ic = ido - i
          val widx1 = i - 1 + iw1
          val widx2 = i - 1 + iw2
          val widx3 = i - 1 + iw3
          w1r = wtable_r(widx1 - 1)
          w1i = wtable_r(widx1)
          w2r = wtable_r(widx2 - 1)
          w2i = wtable_r(widx2)
          w3r = wtable_r(widx3 - 1)
          w3i = wtable_r(widx3)
          val idx9 = in_off + i
          val idx10 = out_off + i
          val idx11 = out_off + ic
          val iidx1 = idx9 + idx1
          val iidx2 = idx9 + idx2
          val iidx3 = idx9 + idx3
          val iidx4 = idx9 + idx4
          val i1i = in(iidx1 - 1)
          val i1r = in(iidx1)
          val i2i = in(iidx2 - 1)
          val i2r = in(iidx2)
          val i3i = in(iidx3 - 1)
          val i3r = in(iidx3)
          val i4i = in(iidx4 - 1)
          val i4r = in(iidx4)
          cr2 = w1r * i2i + w1i * i2r
          ci2 = w1r * i2r - w1i * i2i
          cr3 = w2r * i3i + w2i * i3r
          ci3 = w2r * i3r - w2i * i3i
          cr4 = w3r * i4i + w3i * i4r
          ci4 = w3r * i4r - w3i * i4i
          tr1 = cr2 + cr4
          tr4 = cr4 - cr2
          ti1 = ci2 + ci4
          ti4 = ci2 - ci4
          ti2 = i1r + ci3
          ti3 = i1r - ci3
          tr2 = i1i + cr3
          tr3 = i1i - cr3
          val oidx1 = idx10 + idx5
          val oidx2 = idx11 + idx6
          val oidx3 = idx10 + idx7
          val oidx4 = idx11 + idx8
          out(oidx1 - 1) = tr1 + tr2
          out(oidx4 - 1) = tr2 - tr1
          out(oidx1) = ti1 + ti2
          out(oidx4) = ti1 - ti2
          out(oidx3 - 1) = ti4 + tr3
          out(oidx2 - 1) = tr3 - ti4
          out(oidx3) = tr4 + ti3
          out(oidx2) = tr4 - ti3

          i += 2
        }
      }
      if (ido % 2 == 1) return
    }
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 4 * idx1
      val idx3 = idx1 + idx0
      val idx4 = idx3 + idx0
      val idx5 = idx4 + idx0
      val idx6 = idx2 + ido
      val idx7 = idx6 + ido
      val idx8 = idx7 + ido
      val idx9 = in_off + ido
      val idx10 = out_off + ido
      val i1i = in(idx9 - 1 + idx1)
      val i2i = in(idx9 - 1 + idx3)
      val i3i = in(idx9 - 1 + idx4)
      val i4i = in(idx9 - 1 + idx5)
      ti1 = -hsqt2 * (i2i + i4i)
      tr1 = hsqt2 * (i2i - i4i)
      out(idx10 - 1 + idx2) = tr1 + i1i
      out(idx10 - 1 + idx7) = i1i - tr1
      out(out_off + idx6) = ti1 - i3i
      out(out_off + idx8) = ti1 + i3i
    }
  }

  /*-------------------------------------------------
   radb4: Real FFT's backward processing of factor 4
   -------------------------------------------------*/
  private[fft] def radb4(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    val sqrt2 = 1.41421356237309514547462185873882845
    var i = 0
    var ic = 0
    var ci2 = .0
    var ci3 = .0
    var ci4 = .0
    var cr2 = .0
    var cr3 = .0
    var cr4 = .0
    var ti1 = .0
    var ti2 = .0
    var ti3 = .0
    var ti4 = .0
    var tr1 = .0
    var tr2 = .0
    var tr3 = .0
    var tr4 = .0
    var w1r = .0
    var w1i = .0
    var w2r = .0
    var w2i = .0
    var w3r = .0
    var w3i = .0
    var iw1 = 0
    var iw2 = 0
    var iw3 = 0
    iw1 = offset
    iw2 = iw1 + ido
    iw3 = iw2 + ido
    val idx0 = l1 * ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 4 * idx1
      val idx3 = idx1 + idx0
      val idx4 = idx3 + idx0
      val idx5 = idx4 + idx0
      val idx6 = idx2 + ido
      val idx7 = idx6 + ido
      val idx8 = idx7 + ido
      val i1r = in(in_off + idx2)
      val i2r = in(in_off + idx7)
      val i3r = in(in_off + ido - 1 + idx8)
      val i4r = in(in_off + ido - 1 + idx6)
      tr1 = i1r - i3r
      tr2 = i1r + i3r
      tr3 = i4r + i4r
      tr4 = i2r + i2r
      out(out_off + idx1) = tr2 + tr3
      out(out_off + idx3) = tr1 - tr4
      out(out_off + idx4) = tr2 - tr3
      out(out_off + idx5) = tr1 + tr4
    }
    if (ido < 2) return
    if (ido != 2) {
      for (k <- 0 until l1) {
        val idx1 = k * ido
        val idx2 = idx1 + idx0
        val idx3 = idx2 + idx0
        val idx4 = idx3 + idx0
        val idx5 = 4 * idx1
        val idx6 = idx5 + ido
        val idx7 = idx6 + ido
        val idx8 = idx7 + ido
        i = 2
        while ( {
          i < ido
        }) {
          ic = ido - i
          val widx1 = i - 1 + iw1
          val widx2 = i - 1 + iw2
          val widx3 = i - 1 + iw3
          w1r = wtable_r(widx1 - 1)
          w1i = wtable_r(widx1)
          w2r = wtable_r(widx2 - 1)
          w2i = wtable_r(widx2)
          w3r = wtable_r(widx3 - 1)
          w3i = wtable_r(widx3)
          val idx12 = in_off + i
          val idx13 = in_off + ic
          val idx14 = out_off + i
          val iidx1 = idx12 + idx5
          val iidx2 = idx13 + idx6
          val iidx3 = idx12 + idx7
          val iidx4 = idx13 + idx8
          val i1i = in(iidx1 - 1)
          val i1r = in(iidx1)
          val i2i = in(iidx2 - 1)
          val i2r = in(iidx2)
          val i3i = in(iidx3 - 1)
          val i3r = in(iidx3)
          val i4i = in(iidx4 - 1)
          val i4r = in(iidx4)
          ti1 = i1r + i4r
          ti2 = i1r - i4r
          ti3 = i3r - i2r
          tr4 = i3r + i2r
          tr1 = i1i - i4i
          tr2 = i1i + i4i
          ti4 = i3i - i2i
          tr3 = i3i + i2i
          cr3 = tr2 - tr3
          ci3 = ti2 - ti3
          cr2 = tr1 - tr4
          cr4 = tr1 + tr4
          ci2 = ti1 + ti4
          ci4 = ti1 - ti4
          val oidx1 = idx14 + idx1
          val oidx2 = idx14 + idx2
          val oidx3 = idx14 + idx3
          val oidx4 = idx14 + idx4
          out(oidx1 - 1) = tr2 + tr3
          out(oidx1) = ti2 + ti3
          out(oidx2 - 1) = w1r * cr2 - w1i * ci2
          out(oidx2) = w1r * ci2 + w1i * cr2
          out(oidx3 - 1) = w2r * cr3 - w2i * ci3
          out(oidx3) = w2r * ci3 + w2i * cr3
          out(oidx4 - 1) = w3r * cr4 - w3i * ci4
          out(oidx4) = w3r * ci4 + w3i * cr4

          i += 2
        }
      }
      if (ido % 2 == 1) return
    }
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 4 * idx1
      val idx3 = idx1 + idx0
      val idx4 = idx3 + idx0
      val idx5 = idx4 + idx0
      val idx6 = idx2 + ido
      val idx7 = idx6 + ido
      val idx8 = idx7 + ido
      val idx9 = in_off + ido
      val idx10 = out_off + ido
      val i1r = in(idx9 - 1 + idx2)
      val i2r = in(idx9 - 1 + idx7)
      val i3r = in(in_off + idx6)
      val i4r = in(in_off + idx8)
      ti1 = i3r + i4r
      ti2 = i4r - i3r
      tr1 = i1r - i2r
      tr2 = i1r + i2r
      out(idx10 - 1 + idx1) = tr2 + tr2
      out(idx10 - 1 + idx3) = sqrt2 * (tr1 - ti1)
      out(idx10 - 1 + idx4) = ti2 + ti2
      out(idx10 - 1 + idx5) = -sqrt2 * (tr1 + ti1)
    }
  }

  /*-------------------------------------------------
   radf5: Real FFT's forward processing of factor 5
   -------------------------------------------------*/
  private[fft] def radf5(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    val tr11 = 0.309016994374947451262869435595348477
    val ti11 = 0.951056516295153531181938433292089030
    val tr12 = -0.809016994374947340240566973079694435
    val ti12 = 0.587785252292473248125759255344746634
    var i = 0
    var ic = 0
    var ci2 = .0
    var di2 = .0
    var ci4 = .0
    var ci5 = .0
    var di3 = .0
    var di4 = .0
    var di5 = .0
    var ci3 = .0
    var cr2 = .0
    var cr3 = .0
    var dr2 = .0
    var dr3 = .0
    var dr4 = .0
    var dr5 = .0
    var cr5 = .0
    var cr4 = .0
    var ti2 = .0
    var ti3 = .0
    var ti5 = .0
    var ti4 = .0
    var tr2 = .0
    var tr3 = .0
    var tr4 = .0
    var tr5 = .0
    var w1r = .0
    var w1i = .0
    var w2r = .0
    var w2i = .0
    var w3r = .0
    var w3i = .0
    var w4r = .0
    var w4i = .0
    var iw1 = 0
    var iw2 = 0
    var iw3 = 0
    var iw4 = 0
    iw1 = offset
    iw2 = iw1 + ido
    iw3 = iw2 + ido
    iw4 = iw3 + ido
    val idx0 = l1 * ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 5 * idx1
      val idx3 = idx2 + ido
      val idx4 = idx3 + ido
      val idx5 = idx4 + ido
      val idx6 = idx5 + ido
      val idx7 = idx1 + idx0
      val idx8 = idx7 + idx0
      val idx9 = idx8 + idx0
      val idx10 = idx9 + idx0
      val idx11 = out_off + ido - 1
      val i1r = in(in_off + idx1)
      val i2r = in(in_off + idx7)
      val i3r = in(in_off + idx8)
      val i4r = in(in_off + idx9)
      val i5r = in(in_off + idx10)
      cr2 = i5r + i2r
      ci5 = i5r - i2r
      cr3 = i4r + i3r
      ci4 = i4r - i3r
      out(out_off + idx2) = i1r + cr2 + cr3
      out(idx11 + idx3) = i1r + tr11 * cr2 + tr12 * cr3
      out(out_off + idx4) = ti11 * ci5 + ti12 * ci4
      out(idx11 + idx5) = i1r + tr12 * cr2 + tr11 * cr3
      out(out_off + idx6) = ti12 * ci5 - ti11 * ci4
    }
    if (ido == 1) return
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 5 * idx1
      val idx3 = idx2 + ido
      val idx4 = idx3 + ido
      val idx5 = idx4 + ido
      val idx6 = idx5 + ido
      val idx7 = idx1 + idx0
      val idx8 = idx7 + idx0
      val idx9 = idx8 + idx0
      val idx10 = idx9 + idx0
      i = 2
      while ( {
        i < ido
      }) {
        val widx1 = i - 1 + iw1
        val widx2 = i - 1 + iw2
        val widx3 = i - 1 + iw3
        val widx4 = i - 1 + iw4
        w1r = wtable_r(widx1 - 1)
        w1i = wtable_r(widx1)
        w2r = wtable_r(widx2 - 1)
        w2i = wtable_r(widx2)
        w3r = wtable_r(widx3 - 1)
        w3i = wtable_r(widx3)
        w4r = wtable_r(widx4 - 1)
        w4i = wtable_r(widx4)
        ic = ido - i
        val idx15 = in_off + i
        val idx16 = out_off + i
        val idx17 = out_off + ic
        val iidx1 = idx15 + idx1
        val iidx2 = idx15 + idx7
        val iidx3 = idx15 + idx8
        val iidx4 = idx15 + idx9
        val iidx5 = idx15 + idx10
        val i1i = in(iidx1 - 1)
        val i1r = in(iidx1)
        val i2i = in(iidx2 - 1)
        val i2r = in(iidx2)
        val i3i = in(iidx3 - 1)
        val i3r = in(iidx3)
        val i4i = in(iidx4 - 1)
        val i4r = in(iidx4)
        val i5i = in(iidx5 - 1)
        val i5r = in(iidx5)
        dr2 = w1r * i2i + w1i * i2r
        di2 = w1r * i2r - w1i * i2i
        dr3 = w2r * i3i + w2i * i3r
        di3 = w2r * i3r - w2i * i3i
        dr4 = w3r * i4i + w3i * i4r
        di4 = w3r * i4r - w3i * i4i
        dr5 = w4r * i5i + w4i * i5r
        di5 = w4r * i5r - w4i * i5i
        cr2 = dr2 + dr5
        ci5 = dr5 - dr2
        cr5 = di2 - di5
        ci2 = di2 + di5
        cr3 = dr3 + dr4
        ci4 = dr4 - dr3
        cr4 = di3 - di4
        ci3 = di3 + di4
        tr2 = i1i + tr11 * cr2 + tr12 * cr3
        ti2 = i1r + tr11 * ci2 + tr12 * ci3
        tr3 = i1i + tr12 * cr2 + tr11 * cr3
        ti3 = i1r + tr12 * ci2 + tr11 * ci3
        tr5 = ti11 * cr5 + ti12 * cr4
        ti5 = ti11 * ci5 + ti12 * ci4
        tr4 = ti12 * cr5 - ti11 * cr4
        ti4 = ti12 * ci5 - ti11 * ci4
        val oidx1 = idx16 + idx2
        val oidx2 = idx17 + idx3
        val oidx3 = idx16 + idx4
        val oidx4 = idx17 + idx5
        val oidx5 = idx16 + idx6
        out(oidx1 - 1) = i1i + cr2 + cr3
        out(oidx1) = i1r + ci2 + ci3
        out(oidx3 - 1) = tr2 + tr5
        out(oidx2 - 1) = tr2 - tr5
        out(oidx3) = ti2 + ti5
        out(oidx2) = ti5 - ti2
        out(oidx5 - 1) = tr3 + tr4
        out(oidx4 - 1) = tr3 - tr4
        out(oidx5) = ti3 + ti4
        out(oidx4) = ti4 - ti3

        i += 2
      }
    }
  }

  /*-------------------------------------------------
   radb5: Real FFT's backward processing of factor 5
   -------------------------------------------------*/
  private[fft] def radb5(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    val tr11 = 0.309016994374947451262869435595348477
    val ti11 = 0.951056516295153531181938433292089030
    val tr12 = -0.809016994374947340240566973079694435
    val ti12 = 0.587785252292473248125759255344746634
    var i = 0
    var ic = 0
    var ci2 = .0
    var ci3 = .0
    var ci4 = .0
    var ci5 = .0
    var di3 = .0
    var di4 = .0
    var di5 = .0
    var di2 = .0
    var cr2 = .0
    var cr3 = .0
    var cr5 = .0
    var cr4 = .0
    var ti2 = .0
    var ti3 = .0
    var ti4 = .0
    var ti5 = .0
    var dr3 = .0
    var dr4 = .0
    var dr5 = .0
    var dr2 = .0
    var tr2 = .0
    var tr3 = .0
    var tr4 = .0
    var tr5 = .0
    var w1r = .0
    var w1i = .0
    var w2r = .0
    var w2i = .0
    var w3r = .0
    var w3i = .0
    var w4r = .0
    var w4i = .0
    var iw1 = 0
    var iw2 = 0
    var iw3 = 0
    var iw4 = 0
    iw1 = offset
    iw2 = iw1 + ido
    iw3 = iw2 + ido
    iw4 = iw3 + ido
    val idx0 = l1 * ido
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 5 * idx1
      val idx3 = idx2 + ido
      val idx4 = idx3 + ido
      val idx5 = idx4 + ido
      val idx6 = idx5 + ido
      val idx7 = idx1 + idx0
      val idx8 = idx7 + idx0
      val idx9 = idx8 + idx0
      val idx10 = idx9 + idx0
      val idx11 = in_off + ido - 1
      val i1r = in(in_off + idx2)
      ti5 = 2 * in(in_off + idx4)
      ti4 = 2 * in(in_off + idx6)
      tr2 = 2 * in(idx11 + idx3)
      tr3 = 2 * in(idx11 + idx5)
      cr2 = i1r + tr11 * tr2 + tr12 * tr3
      cr3 = i1r + tr12 * tr2 + tr11 * tr3
      ci5 = ti11 * ti5 + ti12 * ti4
      ci4 = ti12 * ti5 - ti11 * ti4
      out(out_off + idx1) = i1r + tr2 + tr3
      out(out_off + idx7) = cr2 - ci5
      out(out_off + idx8) = cr3 - ci4
      out(out_off + idx9) = cr3 + ci4
      out(out_off + idx10) = cr2 + ci5
    }
    if (ido == 1) return
    for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = 5 * idx1
      val idx3 = idx2 + ido
      val idx4 = idx3 + ido
      val idx5 = idx4 + ido
      val idx6 = idx5 + ido
      val idx7 = idx1 + idx0
      val idx8 = idx7 + idx0
      val idx9 = idx8 + idx0
      val idx10 = idx9 + idx0
      i = 2
      while ( {
        i < ido
      }) {
        ic = ido - i
        val widx1 = i - 1 + iw1
        val widx2 = i - 1 + iw2
        val widx3 = i - 1 + iw3
        val widx4 = i - 1 + iw4
        w1r = wtable_r(widx1 - 1)
        w1i = wtable_r(widx1)
        w2r = wtable_r(widx2 - 1)
        w2i = wtable_r(widx2)
        w3r = wtable_r(widx3 - 1)
        w3i = wtable_r(widx3)
        w4r = wtable_r(widx4 - 1)
        w4i = wtable_r(widx4)
        val idx15 = in_off + i
        val idx16 = in_off + ic
        val idx17 = out_off + i
        val iidx1 = idx15 + idx2
        val iidx2 = idx16 + idx3
        val iidx3 = idx15 + idx4
        val iidx4 = idx16 + idx5
        val iidx5 = idx15 + idx6
        val i1i = in(iidx1 - 1)
        val i1r = in(iidx1)
        val i2i = in(iidx2 - 1)
        val i2r = in(iidx2)
        val i3i = in(iidx3 - 1)
        val i3r = in(iidx3)
        val i4i = in(iidx4 - 1)
        val i4r = in(iidx4)
        val i5i = in(iidx5 - 1)
        val i5r = in(iidx5)
        ti5 = i3r + i2r
        ti2 = i3r - i2r
        ti4 = i5r + i4r
        ti3 = i5r - i4r
        tr5 = i3i - i2i
        tr2 = i3i + i2i
        tr4 = i5i - i4i
        tr3 = i5i + i4i
        cr2 = i1i + tr11 * tr2 + tr12 * tr3
        ci2 = i1r + tr11 * ti2 + tr12 * ti3
        cr3 = i1i + tr12 * tr2 + tr11 * tr3
        ci3 = i1r + tr12 * ti2 + tr11 * ti3
        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4
        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5
        val oidx1 = idx17 + idx1
        val oidx2 = idx17 + idx7
        val oidx3 = idx17 + idx8
        val oidx4 = idx17 + idx9
        val oidx5 = idx17 + idx10
        out(oidx1 - 1) = i1i + tr2 + tr3
        out(oidx1) = i1r + ti2 + ti3
        out(oidx2 - 1) = w1r * dr2 - w1i * di2
        out(oidx2) = w1r * di2 + w1i * dr2
        out(oidx3 - 1) = w2r * dr3 - w2i * di3
        out(oidx3) = w2r * di3 + w2i * dr3
        out(oidx4 - 1) = w3r * dr4 - w3i * di4
        out(oidx4) = w3r * di4 + w3i * dr4
        out(oidx5 - 1) = w4r * dr5 - w4i * di5
        out(oidx5) = w4r * di5 + w4i * dr5

        i += 2
      }
    }
  }

  /*---------------------------------------------------------
   radfg: Real FFT's forward processing of general factor
   --------------------------------------------------------*/
  private[fft] def radfg(ido: Int, ip: Int, l1: Int, idl1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    var idij  = 0
    var j2    = 0
    var ic    = 0
    var jc    = 0
    var lc    = 0
    var is    = 0
    var dc2   = 0.0
    var ai1   = 0.0
    var ai2   = 0.0
    var ar1   = 0.0
    var ar2   = 0.0
    var ds2   = 0.0
    var ar1h  = 0.0
    var ar2h  = 0.0
    var w1r   = 0.0
    var w1i   = 0.0

    val arg   = TWO_PI / ip.toDouble
    val dcp   = cos(arg)
    val dsp   = sin(arg)
    val iw1   = offset
    val ipph  = (ip   + 1) / 2
    val nbd   = (ido  - 1) / 2

    if (ido != 1) {
      for (ik <- 0 until idl1) {
        out(out_off + ik) = in(in_off + ik)
      }
      for (j <- 1 until ip) {
        val idx1 = j * l1 * ido
        for (k <- 0 until l1) {
          val idx2 = k * ido + idx1
          out(out_off + idx2) = in(in_off + idx2)
        }
      }
      if (nbd <= l1) {
        is = -ido
        for (j <- 1 until ip) {
          is += ido
          idij = is - 1
          val idx1 = j * l1 * ido
          var i = 2
          while (i < ido) {
            idij += 2
            val idx2 = idij + iw1
            val idx4 = in_off + i
            val idx5 = out_off + i
            w1r = wtable_r(idx2 - 1)
            w1i = wtable_r(idx2)
            for (k <- 0 until l1) {
              val idx3  = k * ido + idx1
              val oidx1 = idx5 + idx3
              val iidx1 = idx4 + idx3
              val i1i   = in(iidx1 - 1)
              val i1r   = in(iidx1)
              out(oidx1 - 1 ) = w1r * i1i + w1i * i1r
              out(oidx1     ) = w1r * i1r - w1i * i1i
            }

            i += 2
          }
        }
      }
      else {
        is = -ido
        for (j <- 1 until ip) {
          is += ido
          val idx1 = j * l1 * ido
          for (k <- 0 until l1) {
            idij = is - 1
            val idx3 = k * ido + idx1
            var i = 2
            while (i < ido) {
              idij += 2
              val idx2  = idij + iw1
              w1r       = wtable_r(idx2 - 1)
              w1i       = wtable_r(idx2)
              val oidx1 = out_off + i + idx3
              val iidx1 = in_off + i + idx3
              val i1i   = in(iidx1 - 1)
              val i1r   = in(iidx1)
              out(oidx1 - 1 ) = w1r * i1i + w1i * i1r
              out(oidx1     ) = w1r * i1r - w1i * i1i

              i += 2
            }
          }
        }
      }
      if (nbd >= l1) for (j <- 1 until ipph) {
        jc = ip - j
        val idx1 = j * l1 * ido
        val idx2 = jc * l1 * ido
        for (k <- 0 until l1) {
          val idx3 = k * ido + idx1
          val idx4 = k * ido + idx2
          var i = 2
          while (i < ido) {
            val idx5  = in_off + i
            val idx6  = out_off + i
            val iidx1 = idx5 + idx3
            val iidx2 = idx5 + idx4
            val oidx1 = idx6 + idx3
            val oidx2 = idx6 + idx4
            val o1i   = out(oidx1 - 1)
            val o1r   = out(oidx1)
            val o2i   = out(oidx2 - 1)
            val o2r   = out(oidx2)
            in(iidx1 - 1) = o1i + o2i
            in(iidx1    ) = o1r + o2r
            in(iidx2 - 1) = o1r - o2r
            in(iidx2    ) = o2i - o1i

            i += 2
          }
        }
      }
      else for (j <- 1 until ipph) {
        jc = ip - j
        val idx1 = j * l1 * ido
        val idx2 = jc * l1 * ido
        var i = 2
        while (i < ido) {
          val idx5 = in_off + i
          val idx6 = out_off + i
          for (k <- 0 until l1) {
            val idx3  = k * ido + idx1
            val idx4  = k * ido + idx2
            val iidx1 = idx5 + idx3
            val iidx2 = idx5 + idx4
            val oidx1 = idx6 + idx3
            val oidx2 = idx6 + idx4
            val o1i   = out(oidx1 - 1)
            val o1r   = out(oidx1)
            val o2i   = out(oidx2 - 1)
            val o2r   = out(oidx2)
            in(iidx1 - 1) = o1i + o2i
            in(iidx1    ) = o1r + o2r
            in(iidx2 - 1) = o1r - o2r
            in(iidx2    ) = o2i - o1i
          }

          i += 2
        }
      }
    }
    else System.arraycopy(out, out_off, in, in_off, idl1)
    for (j <- 1 until ipph) {
      jc = ip - j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      for (k <- 0 until l1) {
        val idx3  = k * ido + idx1
        val idx4  = k * ido + idx2
        val oidx1 = out_off + idx3
        val oidx2 = out_off + idx4
        val o1r   = out(oidx1)
        val o2r   = out(oidx2)
        in(in_off + idx3) = o1r + o2r
        in(in_off + idx4) = o2r - o1r
      }
    }
    ar1 = 1
    ai1 = 0
    val idx0 = (ip - 1) * idl1
    for (l <- 1 until ipph) {
      lc    = ip - l
      ar1h  = dcp * ar1 - dsp * ai1
      ai1   = dcp * ai1 + dsp * ar1
      ar1   = ar1h
      val idx1 = l * idl1
      val idx2 = lc * idl1
      for (ik <- 0 until idl1) {
        val idx3 = out_off + ik
        val idx4 = in_off + ik
        out(idx3 + idx1) = in(idx4) + ar1 * in(idx4 + idl1)
        out(idx3 + idx2) = ai1 * in(idx4 + idx0)
      }
      dc2 = ar1
      ds2 = ai1
      ar2 = ar1
      ai2 = ai1
      for (j <- 2 until ipph) {
        jc    = ip - j
        ar2h  = dc2 * ar2 - ds2 * ai2
        ai2   = dc2 * ai2 + ds2 * ar2
        ar2   = ar2h
        val idx3 = j * idl1
        val idx4 = jc * idl1
        for (ik <- 0 until idl1) {
          val idx5 = out_off + ik
          val idx6 = in_off + ik
          out(idx5 + idx1) += ar2 * in(idx6 + idx3)
          out(idx5 + idx2) += ai2 * in(idx6 + idx4)
        }
      }
    }
    for (j <- 1 until ipph) {
      val idx1 = j * idl1
      for (ik <- 0 until idl1) {
        out(out_off + ik) += in(in_off + ik + idx1)
      }
    }
    if (ido >= l1) for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = idx1 * ip
      for (i <- 0 until ido) {
        in(in_off + i + idx2) = out(out_off + i + idx1)
      }
    }
    else for (i <- 0 until ido) {
      for (k <- 0 until l1) {
        val idx1 = k * ido
        in(in_off + i + idx1 * ip) = out(out_off + i + idx1)
      }
    }
    val idx01 = ip * ido
    for (j <- 1 until ipph) {
      jc = ip - j
      j2 = 2 * j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      val idx3 = j2 * ido
      for (k <- 0 until l1) {
        val idx4 = k * ido
        val idx5 = idx4 + idx1
        val idx6 = idx4 + idx2
        val idx7 = k * idx01
        in(in_off + ido - 1 + idx3 - ido + idx7) = out(out_off + idx5)
        in(in_off + idx3 + idx7) = out(out_off + idx6)
      }
    }
    if (ido == 1) return
    if (nbd >= l1) for (j <- 1 until ipph) {
      jc = ip - j
      j2 = 2 * j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      val idx3 = j2 * ido
      for (k <- 0 until l1) {
        val idx4 = k * idx01
        val idx5 = k * ido
        var i = 2
        while (i < ido) {
          ic = ido - i
          val idx6  = in_off + i
          val idx7  = in_off + ic
          val idx8  = out_off + i
          val iidx1 = idx6 + idx3 + idx4
          val iidx2 = idx7 + idx3 - ido + idx4
          val oidx1 = idx8 + idx5 + idx1
          val oidx2 = idx8 + idx5 + idx2
          val o1i   = out(oidx1 - 1)
          val o1r   = out(oidx1)
          val o2i   = out(oidx2 - 1)
          val o2r   = out(oidx2)
          in(iidx1 - 1) = o1i + o2i
          in(iidx2 - 1) = o1i - o2i
          in(iidx1    ) = o1r + o2r
          in(iidx2    ) = o2r - o1r

          i += 2
        }
      }
    }
    else for (j <- 1 until ipph) {
      jc = ip - j
      j2 = 2 * j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      val idx3 = j2 * ido
      var i = 2
      while (i < ido) {
        ic = ido - i
        val idx6 = in_off + i
        val idx7 = in_off + ic
        val idx8 = out_off + i
        for (k <- 0 until l1) {
          val idx4  = k * idx01
          val idx5  = k * ido
          val iidx1 = idx6 + idx3 + idx4
          val iidx2 = idx7 + idx3 - ido + idx4
          val oidx1 = idx8 + idx5 + idx1
          val oidx2 = idx8 + idx5 + idx2
          val o1i   = out(oidx1 - 1)
          val o1r   = out(oidx1)
          val o2i   = out(oidx2 - 1)
          val o2r   = out(oidx2)
          in(iidx1 - 1) = o1i + o2i
          in(iidx2 - 1) = o1i - o2i
          in(iidx1    ) = o1r + o2r
          in(iidx2    ) = o2r - o1r
        }

        i += 2
      }
    }
  }

  /*---------------------------------------------------------
   radbg: Real FFT's backward processing of general factor
   --------------------------------------------------------*/
  private[fft] def radbg(ido: Int, ip: Int, l1: Int, idl1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int): Unit = {
    val arg   = TWO_PI / ip.toDouble
    val dcp   = cos(arg)
    val dsp   = sin(arg)
    val iw1   = offset
    val nbd   = (ido  - 1) / 2
    val ipph  = (ip   + 1) / 2

    var idij  = 0
    var j2    = 0
    var ic    = 0
    var jc    = 0
    var lc    = 0
    var is    = 0
    var dc2   = 0.0
    var ai1   = 0.0
    var ai2   = 0.0
    var ar1   = 0.0
    var ar2   = 0.0
    var ds2   = 0.0
    var w1r   = 0.0
    var w1i   = 0.0
    var ar1h  = 0.0
    var ar2h  = 0.0

    val idx0 = ip * ido
    if (ido >= l1) for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = k * idx0
      for (i <- 0 until ido) {
        out(out_off + i + idx1) = in(in_off + i + idx2)
      }
    }
    else for (i <- 0 until ido) {
      val idx1 = out_off + i
      val idx2 = in_off + i
      for (k <- 0 until l1) {
        out(idx1 + k * ido) = in(idx2 + k * idx0)
      }
    }
    val iidx0 = in_off + ido - 1
    for (j <- 1 until ipph) {
      jc = ip - j
      j2 = 2 * j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      val idx3 = j2 * ido
      for (k <- 0 until l1) {
        val idx4  = k * ido
        val idx5  = idx4 * ip
        val iidx1 = iidx0 + idx3 + idx5 - ido
        val iidx2 = in_off + idx3 + idx5
        val i1r   = in(iidx1)
        val i2r   = in(iidx2)
        out(out_off + idx4 + idx1) = i1r + i1r
        out(out_off + idx4 + idx2) = i2r + i2r
      }
    }
    if (ido != 1) if (nbd >= l1) for (j <- 1 until ipph) {
      jc = ip - j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      val idx3 = 2 * j * ido
      for (k <- 0 until l1) {
        val idx4 = k * ido + idx1
        val idx5 = k * ido + idx2
        val idx6 = k * ip * ido + idx3
        var i = 2
        while (i < ido) {
          ic = ido - i
          val idx7  = out_off + i
          val idx8  = in_off + ic
          val idx9  = in_off + i
          val oidx1 = idx7 + idx4
          val oidx2 = idx7 + idx5
          val iidx1 = idx9 + idx6
          val iidx2 = idx8 + idx6 - ido
          val a1i   = in(iidx1 - 1)
          val a1r   = in(iidx1)
          val a2i   = in(iidx2 - 1)
          val a2r   = in(iidx2)
          out(oidx1 - 1) = a1i + a2i
          out(oidx2 - 1) = a1i - a2i
          out(oidx1   ) = a1r - a2r
          out(oidx2   ) = a1r + a2r

          i += 2
        }
      }
    }
    else for (j <- 1 until ipph) {
      jc = ip - j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      val idx3 = 2 * j * ido
      var i = 2
      while (i < ido) {
        ic = ido - i
        val idx7 = out_off + i
        val idx8 = in_off + ic
        val idx9 = in_off + i
        for (k <- 0 until l1) {
          val idx4  = k * ido + idx1
          val idx5  = k * ido + idx2
          val idx6  = k * ip * ido + idx3
          val oidx1 = idx7 + idx4
          val oidx2 = idx7 + idx5
          val iidx1 = idx9 + idx6
          val iidx2 = idx8 + idx6 - ido
          val a1i   = in(iidx1 - 1)
          val a1r   = in(iidx1)
          val a2i   = in(iidx2 - 1)
          val a2r   = in(iidx2)
          out(oidx1 - 1 ) = a1i + a2i
          out(oidx2 - 1 ) = a1i - a2i
          out(oidx1     ) = a1r - a2r
          out(oidx2     ) = a1r + a2r
        }

        i += 2
      }
    }
    ar1 = 1
    ai1 = 0
    val idx01 = (ip - 1) * idl1
    for (l <- 1 until ipph) {
      lc    = ip - l
      ar1h  = dcp * ar1 - dsp * ai1
      ai1   = dcp * ai1 + dsp * ar1
      ar1   = ar1h
      val idx1 = l * idl1
      val idx2 = lc * idl1
      for (ik <- 0 until idl1) {
        val idx3 = in_off + ik
        val idx4 = out_off + ik
        in(idx3 + idx1) = out(idx4) + ar1 * out(idx4 + idl1)
        in(idx3 + idx2) = ai1 * out(idx4 + idx01)
      }
      dc2 = ar1
      ds2 = ai1
      ar2 = ar1
      ai2 = ai1
      for (j <- 2 until ipph) {
        jc    = ip - j
        ar2h  = dc2 * ar2 - ds2 * ai2
        ai2   = dc2 * ai2 + ds2 * ar2
        ar2   = ar2h
        val idx5 = j * idl1
        val idx6 = jc * idl1
        for (ik <- 0 until idl1) {
          val idx7 = in_off + ik
          val idx8 = out_off + ik
          in(idx7 + idx1) += ar2 * out(idx8 + idx5)
          in(idx7 + idx2) += ai2 * out(idx8 + idx6)
        }
      }
    }
    for (j <- 1 until ipph) {
      val idx1 = j * idl1
      for (ik <- 0 until idl1) {
        val idx2 = out_off + ik
        out(idx2) += out(idx2 + idx1)
      }
    }
    for (j <- 1 until ipph) {
      jc = ip - j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      for (k <- 0 until l1) {
        val idx3  = k * ido
        val oidx1 = out_off + idx3
        val iidx1 = in_off + idx3 + idx1
        val iidx2 = in_off + idx3 + idx2
        val i1r   = in(iidx1)
        val i2r   = in(iidx2)
        out(oidx1 + idx1) = i1r - i2r
        out(oidx1 + idx2) = i1r + i2r
      }
    }
    if (ido == 1) return
    if (nbd >= l1) for (j <- 1 until ipph) {
      jc = ip - j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      for (k <- 0 until l1) {
        val idx3 = k * ido
        var i = 2
        while (i < ido) {
          val idx4 = out_off + i
          val idx5 = in_off  + i
          val oidx1 = idx4 + idx3 + idx1
          val oidx2 = idx4 + idx3 + idx2
          val iidx1 = idx5 + idx3 + idx1
          val iidx2 = idx5 + idx3 + idx2
          val i1i = in(iidx1 - 1)
          val i1r = in(iidx1)
          val i2i = in(iidx2 - 1)
          val i2r = in(iidx2)
          out(oidx1 - 1 ) = i1i - i2r
          out(oidx2 - 1 ) = i1i + i2r
          out(oidx1     ) = i1r + i2i
          out(oidx2     ) = i1r - i2i

          i += 2
        }
      }
    }
    else for (j <- 1 until ipph) {
      jc = ip - j
      val idx1 = j * l1 * ido
      val idx2 = jc * l1 * ido
      var i = 2
      while (i < ido) {
        val idx4 = out_off + i
        val idx5 = in_off  + i
        for (k <- 0 until l1) {
          val idx3  = k * ido
          val oidx1 = idx4 + idx3 + idx1
          val oidx2 = idx4 + idx3 + idx2
          val iidx1 = idx5 + idx3 + idx1
          val iidx2 = idx5 + idx3 + idx2
          val i1i   = in(iidx1 - 1)
          val i1r   = in(iidx1)
          val i2i   = in(iidx2 - 1)
          val i2r   = in(iidx2)
          out(oidx1 - 1 ) = i1i - i2r
          out(oidx2 - 1 ) = i1i + i2r
          out(oidx1     ) = i1r + i2i
          out(oidx2     ) = i1r - i2i
        }

        i += 2
      }
    }
    System.arraycopy(out, out_off, in, in_off, idl1)
    for (j <- 1 until ip) {
      val idx1 = j * l1 * ido
      for (k <- 0 until l1) {
        val idx2 = k * ido + idx1
        in(in_off + idx2) = out(out_off + idx2)
      }
    }
    if (nbd <= l1) {
      is = -ido
      for (j <- 1 until ip) {
        is += ido
        idij = is - 1
        val idx1 = j * l1 * ido
        var i = 2
        while (i < ido) {
          idij += 2
          val idx2 = idij + iw1
          w1r = wtable_r(idx2 - 1)
          w1i = wtable_r(idx2)
          val idx4 = in_off + i
          val idx5 = out_off + i
          for (k <- 0 until l1) {
            val idx3  = k * ido + idx1
            val iidx1 = idx4 + idx3
            val oidx1 = idx5 + idx3
            val o1i   = out(oidx1 - 1)
            val o1r   = out(oidx1)
            in(iidx1 - 1) = w1r * o1i - w1i * o1r
            in(iidx1    ) = w1r * o1r + w1i * o1i
          }

          i += 2
        }
      }
    }
    else {
      is = -ido
      for (j <- 1 until ip) {
        is += ido
        val idx1 = j * l1 * ido
        for (k <- 0 until l1) {
          idij = is - 1
          val idx3 = k * ido + idx1
          var i = 2
          while (i < ido) {
            idij += 2
            val idx2 = idij + iw1
            w1r = wtable_r(idx2 - 1)
            w1i = wtable_r(idx2)
            val idx4  = in_off + i
            val idx5  = out_off + i
            val iidx1 = idx4 + idx3
            val oidx1 = idx5 + idx3
            val o1i   = out(oidx1 - 1)
            val o1r   = out(oidx1)
            in(iidx1 - 1) = w1r * o1i - w1i * o1r
            in(iidx1    ) = w1r * o1r + w1i * o1i

            i += 2
          }
        }
      }
    }
  }

  /*---------------------------------------------------------
     cfftf1: further processing of Complex forward FFT
     --------------------------------------------------------*/
  private[fft] def cfftf(a: Array[Double], offa: Int, isign: Int): Unit = {
    var idot = 0
    var l1 = 0
    var l2 = 0
    var na = 0
    var nf = 0
    var ipll = 0
    var iw = 0
    var ido = 0
    var idl1 = 0
    val nac = new Array[Int](1)
    val twon = 2 * n
    var iw1 = 0
    var iw2 = 0
    val ch = new Array[Double](twon)
    iw1 = twon
    iw2 = 4 * n
    nac(0) = 0
    nf = wtable(1 + iw2).toInt
    na = 0
    l1 = 1
    iw = iw1
    for (k1 <- 2 to nf + 1) {
      ipll = wtable(k1 + iw2).toInt
      l2 = ipll * l1
      ido = n / l2
      idot = ido + ido
      idl1 = idot * l1
      ipll match {
        case 4 =>
          if (na == 0) passf4(idot, l1, a, offa, ch, 0, iw, isign)
          else passf4(idot, l1, ch, 0, a, offa, iw, isign)
          na = 1 - na

        case 2 =>
          if (na == 0) passf2(idot, l1, a, offa, ch, 0, iw, isign)
          else passf2(idot, l1, ch, 0, a, offa, iw, isign)
          na = 1 - na

        case 3 =>
          if (na == 0) passf3(idot, l1, a, offa, ch, 0, iw, isign)
          else passf3(idot, l1, ch, 0, a, offa, iw, isign)
          na = 1 - na

        case 5 =>
          if (na == 0) passf5(idot, l1, a, offa, ch, 0, iw, isign)
          else passf5(idot, l1, ch, 0, a, offa, iw, isign)
          na = 1 - na

        case _ =>
          if (na == 0) passfg(nac, idot, ipll, l1, idl1, a, offa, ch, 0, iw, isign)
          else passfg(nac, idot, ipll, l1, idl1, ch, 0, a, offa, iw, isign)
          if (nac(0) != 0) na = 1 - na

      }
      l1 = l2
      iw += (ipll - 1) * idot
    }
    if (na == 0) return
    System.arraycopy(ch, 0, a, offa, twon)
  }

  /*----------------------------------------------------------------------
   passf2: Complex FFT's forward/backward processing of factor 2;
   isign is +1 for backward and -1 for forward transforms
   ----------------------------------------------------------------------*/
  private[fft] def passf2(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int, isign: Int): Unit = {
    var t1i = .0
    var t1r = .0
    var iw1 = 0
    iw1 = offset
    val idx = ido * l1
    if (ido <= 2) for (k <- 0 until l1) {
      val idx0 = k * ido
      val iidx1 = in_off + 2 * idx0
      val iidx2 = iidx1 + ido
      val a1r = in(iidx1)
      val a1i = in(iidx1 + 1)
      val a2r = in(iidx2)
      val a2i = in(iidx2 + 1)
      val oidx1 = out_off + idx0
      val oidx2 = oidx1 + idx
      out(oidx1) = a1r + a2r
      out(oidx1 + 1) = a1i + a2i
      out(oidx2) = a1r - a2r
      out(oidx2 + 1) = a1i - a2i
    }
    else for (k <- 0 until l1) {
      var i = 0
      while ( {
        i < ido - 1
      }) {
        val idx0 = k * ido
        val iidx1 = in_off + i + 2 * idx0
        val iidx2 = iidx1 + ido
        val i1r = in(iidx1)
        val i1i = in(iidx1 + 1)
        val i2r = in(iidx2)
        val i2i = in(iidx2 + 1)
        val widx1 = i + iw1
        val w1r = wtable(widx1)
        val w1i = isign * wtable(widx1 + 1)
        t1r = i1r - i2r
        t1i = i1i - i2i
        val oidx1 = out_off + i + idx0
        val oidx2 = oidx1 + idx
        out(oidx1) = i1r + i2r
        out(oidx1 + 1) = i1i + i2i
        out(oidx2) = w1r * t1r - w1i * t1i
        out(oidx2 + 1) = w1r * t1i + w1i * t1r

        i += 2
      }
    }
  }

  /*----------------------------------------------------------------------
   passf3: Complex FFT's forward/backward processing of factor 3;
   isign is +1 for backward and -1 for forward transforms
   ----------------------------------------------------------------------*/
  private[fft] def passf3(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int, isign: Int): Unit = {
    val taur = -0.5
    val taui = 0.866025403784438707610604524234076962
    var ci2 = .0
    var ci3 = .0
    var di2 = .0
    var di3 = .0
    var cr2 = .0
    var cr3 = .0
    var dr2 = .0
    var dr3 = .0
    var ti2 = .0
    var tr2 = .0
    var iw1 = 0
    var iw2 = 0
    iw1 = offset
    iw2 = iw1 + ido
    val idxt = l1 * ido
    if (ido == 2) for (k <- 1 to l1) {
      val iidx1 = in_off + (3 * k - 2) * ido
      val iidx2 = iidx1 + ido
      val iidx3 = iidx1 - ido
      val i1r = in(iidx1)
      val i1i = in(iidx1 + 1)
      val i2r = in(iidx2)
      val i2i = in(iidx2 + 1)
      val i3r = in(iidx3)
      val i3i = in(iidx3 + 1)
      tr2 = i1r + i2r
      cr2 = i3r + taur * tr2
      ti2 = i1i + i2i
      ci2 = i3i + taur * ti2
      cr3 = isign * taui * (i1r - i2r)
      ci3 = isign * taui * (i1i - i2i)
      val oidx1 = out_off + (k - 1) * ido
      val oidx2 = oidx1 + idxt
      val oidx3 = oidx2 + idxt
      out(oidx1) = in(iidx3) + tr2
      out(oidx1 + 1) = i3i + ti2
      out(oidx2) = cr2 - ci3
      out(oidx2 + 1) = ci2 + cr3
      out(oidx3) = cr2 + ci3
      out(oidx3 + 1) = ci2 - cr3
    }
    else for (k <- 1 to l1) {
      val idx1 = in_off + (3 * k - 2) * ido
      val idx2 = out_off + (k - 1) * ido
      var i = 0
      while ( {
        i < ido - 1
      }) {
        val iidx1 = i + idx1
        val iidx2 = iidx1 + ido
        val iidx3 = iidx1 - ido
        val a1r = in(iidx1)
        val a1i = in(iidx1 + 1)
        val a2r = in(iidx2)
        val a2i = in(iidx2 + 1)
        val a3r = in(iidx3)
        val a3i = in(iidx3 + 1)
        tr2 = a1r + a2r
        cr2 = a3r + taur * tr2
        ti2 = a1i + a2i
        ci2 = a3i + taur * ti2
        cr3 = isign * taui * (a1r - a2r)
        ci3 = isign * taui * (a1i - a2i)
        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3
        val widx1 = i + iw1
        val widx2 = i + iw2
        val w1r = wtable(widx1)
        val w1i = isign * wtable(widx1 + 1)
        val w2r = wtable(widx2)
        val w2i = isign * wtable(widx2 + 1)
        val oidx1 = i + idx2
        val oidx2 = oidx1 + idxt
        val oidx3 = oidx2 + idxt
        out(oidx1) = a3r + tr2
        out(oidx1 + 1) = a3i + ti2
        out(oidx2) = w1r * dr2 - w1i * di2
        out(oidx2 + 1) = w1r * di2 + w1i * dr2
        out(oidx3) = w2r * dr3 - w2i * di3
        out(oidx3 + 1) = w2r * di3 + w2i * dr3

        i += 2
      }
    }
  }

  /*----------------------------------------------------------------------
   passf4: Complex FFT's forward/backward processing of factor 4;
   isign is +1 for backward and -1 for forward transforms
   ----------------------------------------------------------------------*/
  private[fft] def passf4(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int, isign: Int): Unit = {
    var ci2 = .0
    var ci3 = .0
    var ci4 = .0
    var cr2 = .0
    var cr3 = .0
    var cr4 = .0
    var ti1 = .0
    var ti2 = .0
    var ti3 = .0
    var ti4 = .0
    var tr1 = .0
    var tr2 = .0
    var tr3 = .0
    var tr4 = .0
    var iw1 = 0
    var iw2 = 0
    var iw3 = 0
    iw1 = offset
    iw2 = iw1 + ido
    iw3 = iw2 + ido
    val idx0 = l1 * ido
    if (ido == 2) for (k <- 0 until l1) {
      val idxt1 = k * ido
      val iidx1 = in_off + 4 * idxt1 + 1
      val iidx2 = iidx1 + ido
      val iidx3 = iidx2 + ido
      val iidx4 = iidx3 + ido
      val i1i = in(iidx1 - 1)
      val i1r = in(iidx1)
      val i2i = in(iidx2 - 1)
      val i2r = in(iidx2)
      val i3i = in(iidx3 - 1)
      val i3r = in(iidx3)
      val i4i = in(iidx4 - 1)
      val i4r = in(iidx4)
      ti1 = i1r - i3r
      ti2 = i1r + i3r
      tr4 = i4r - i2r
      ti3 = i2r + i4r
      tr1 = i1i - i3i
      tr2 = i1i + i3i
      ti4 = i2i - i4i
      tr3 = i2i + i4i
      val oidx1 = out_off + idxt1
      val oidx2 = oidx1 + idx0
      val oidx3 = oidx2 + idx0
      val oidx4 = oidx3 + idx0
      out(oidx1) = tr2 + tr3
      out(oidx1 + 1) = ti2 + ti3
      out(oidx2) = tr1 + isign * tr4
      out(oidx2 + 1) = ti1 + isign * ti4
      out(oidx3) = tr2 - tr3
      out(oidx3 + 1) = ti2 - ti3
      out(oidx4) = tr1 - isign * tr4
      out(oidx4 + 1) = ti1 - isign * ti4
    }
    else for (k <- 0 until l1) {
      val idx1 = k * ido
      val idx2 = in_off + 1 + 4 * idx1
      var i = 0
      while ( {
        i < ido - 1
      }) {
        val iidx1 = i + idx2
        val iidx2 = iidx1 + ido
        val iidx3 = iidx2 + ido
        val iidx4 = iidx3 + ido
        val i1i = in(iidx1 - 1)
        val i1r = in(iidx1)
        val i2i = in(iidx2 - 1)
        val i2r = in(iidx2)
        val i3i = in(iidx3 - 1)
        val i3r = in(iidx3)
        val i4i = in(iidx4 - 1)
        val i4r = in(iidx4)
        ti1 = i1r - i3r
        ti2 = i1r + i3r
        ti3 = i2r + i4r
        tr4 = i4r - i2r
        tr1 = i1i - i3i
        tr2 = i1i + i3i
        ti4 = i2i - i4i
        tr3 = i2i + i4i
        cr3 = tr2 - tr3
        ci3 = ti2 - ti3
        cr2 = tr1 + isign * tr4
        cr4 = tr1 - isign * tr4
        ci2 = ti1 + isign * ti4
        ci4 = ti1 - isign * ti4
        val widx1 = i + iw1
        val widx2 = i + iw2
        val widx3 = i + iw3
        val w1r = wtable(widx1)
        val w1i = isign * wtable(widx1 + 1)
        val w2r = wtable(widx2)
        val w2i = isign * wtable(widx2 + 1)
        val w3r = wtable(widx3)
        val w3i = isign * wtable(widx3 + 1)
        val oidx1 = out_off + i + idx1
        val oidx2 = oidx1 + idx0
        val oidx3 = oidx2 + idx0
        val oidx4 = oidx3 + idx0
        out(oidx1) = tr2 + tr3
        out(oidx1 + 1) = ti2 + ti3
        out(oidx2) = w1r * cr2 - w1i * ci2
        out(oidx2 + 1) = w1r * ci2 + w1i * cr2
        out(oidx3) = w2r * cr3 - w2i * ci3
        out(oidx3 + 1) = w2r * ci3 + w2i * cr3
        out(oidx4) = w3r * cr4 - w3i * ci4
        out(oidx4 + 1) = w3r * ci4 + w3i * cr4

        i += 2
      }
    }
  }

  /*----------------------------------------------------------------------
   passf5: Complex FFT's forward/backward processing of factor 5;
   isign is +1 for backward and -1 for forward transforms
   ----------------------------------------------------------------------*/
  private[fft] def passf5(ido: Int, l1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int, isign: Int): Unit = {
    val tr11 = 0.309016994374947451262869435595348477
    val ti11 = 0.951056516295153531181938433292089030
    val tr12 = -0.809016994374947340240566973079694435
    val ti12 = 0.587785252292473248125759255344746634
    var ci2 = .0
    var ci3 = .0
    var ci4 = .0
    var ci5 = .0
    var di3 = .0
    var di4 = .0
    var di5 = .0
    var di2 = .0
    var cr2 = .0
    var cr3 = .0
    var cr5 = .0
    var cr4 = .0
    var ti2 = .0
    var ti3 = .0
    var ti4 = .0
    var ti5 = .0
    var dr3 = .0
    var dr4 = .0
    var dr5 = .0
    var dr2 = .0
    var tr2 = .0
    var tr3 = .0
    var tr4 = .0
    var tr5 = .0
    var iw1 = 0
    var iw2 = 0
    var iw3 = 0
    var iw4 = 0
    iw1 = offset
    iw2 = iw1 + ido
    iw3 = iw2 + ido
    iw4 = iw3 + ido
    val idx0 = l1 * ido
    if (ido == 2) for (k <- 1 to l1) {
      val iidx1 = in_off + (5 * k - 4) * ido + 1
      val iidx2 = iidx1 + ido
      val iidx3 = iidx1 - ido
      val iidx4 = iidx2 + ido
      val iidx5 = iidx4 + ido
      val i1i = in(iidx1 - 1)
      val i1r = in(iidx1)
      val i2i = in(iidx2 - 1)
      val i2r = in(iidx2)
      val i3i = in(iidx3 - 1)
      val i3r = in(iidx3)
      val i4i = in(iidx4 - 1)
      val i4r = in(iidx4)
      val i5i = in(iidx5 - 1)
      val i5r = in(iidx5)
      ti5 = i1r - i5r
      ti2 = i1r + i5r
      ti4 = i2r - i4r
      ti3 = i2r + i4r
      tr5 = i1i - i5i
      tr2 = i1i + i5i
      tr4 = i2i - i4i
      tr3 = i2i + i4i
      cr2 = i3i + tr11 * tr2 + tr12 * tr3
      ci2 = i3r + tr11 * ti2 + tr12 * ti3
      cr3 = i3i + tr12 * tr2 + tr11 * tr3
      ci3 = i3r + tr12 * ti2 + tr11 * ti3
      cr5 = isign * (ti11 * tr5 + ti12 * tr4)
      ci5 = isign * (ti11 * ti5 + ti12 * ti4)
      cr4 = isign * (ti12 * tr5 - ti11 * tr4)
      ci4 = isign * (ti12 * ti5 - ti11 * ti4)
      val oidx1 = out_off + (k - 1) * ido
      val oidx2 = oidx1 + idx0
      val oidx3 = oidx2 + idx0
      val oidx4 = oidx3 + idx0
      val oidx5 = oidx4 + idx0
      out(oidx1) = i3i + tr2 + tr3
      out(oidx1 + 1) = i3r + ti2 + ti3
      out(oidx2) = cr2 - ci5
      out(oidx2 + 1) = ci2 + cr5
      out(oidx3) = cr3 - ci4
      out(oidx3 + 1) = ci3 + cr4
      out(oidx4) = cr3 + ci4
      out(oidx4 + 1) = ci3 - cr4
      out(oidx5) = cr2 + ci5
      out(oidx5 + 1) = ci2 - cr5
    }
    else for (k <- 1 to l1) {
      val idx1 = in_off + 1 + (k * 5 - 4) * ido
      val idx2 = out_off + (k - 1) * ido
      var i = 0
      while ( {
        i < ido - 1
      }) {
        val iidx1 = i + idx1
        val iidx2 = iidx1 + ido
        val iidx3 = iidx1 - ido
        val iidx4 = iidx2 + ido
        val iidx5 = iidx4 + ido
        val i1i = in(iidx1 - 1)
        val i1r = in(iidx1)
        val i2i = in(iidx2 - 1)
        val i2r = in(iidx2)
        val i3i = in(iidx3 - 1)
        val i3r = in(iidx3)
        val i4i = in(iidx4 - 1)
        val i4r = in(iidx4)
        val i5i = in(iidx5 - 1)
        val i5r = in(iidx5)
        ti5 = i1r - i5r
        ti2 = i1r + i5r
        ti4 = i2r - i4r
        ti3 = i2r + i4r
        tr5 = i1i - i5i
        tr2 = i1i + i5i
        tr4 = i2i - i4i
        tr3 = i2i + i4i
        cr2 = i3i + tr11 * tr2 + tr12 * tr3
        ci2 = i3r + tr11 * ti2 + tr12 * ti3
        cr3 = i3i + tr12 * tr2 + tr11 * tr3
        ci3 = i3r + tr12 * ti2 + tr11 * ti3
        cr5 = isign * (ti11 * tr5 + ti12 * tr4)
        ci5 = isign * (ti11 * ti5 + ti12 * ti4)
        cr4 = isign * (ti12 * tr5 - ti11 * tr4)
        ci4 = isign * (ti12 * ti5 - ti11 * ti4)
        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5
        val widx1 = i + iw1
        val widx2 = i + iw2
        val widx3 = i + iw3
        val widx4 = i + iw4
        val w1r = wtable(widx1)
        val w1i = isign * wtable(widx1 + 1)
        val w2r = wtable(widx2)
        val w2i = isign * wtable(widx2 + 1)
        val w3r = wtable(widx3)
        val w3i = isign * wtable(widx3 + 1)
        val w4r = wtable(widx4)
        val w4i = isign * wtable(widx4 + 1)
        val oidx1 = i + idx2
        val oidx2 = oidx1 + idx0
        val oidx3 = oidx2 + idx0
        val oidx4 = oidx3 + idx0
        val oidx5 = oidx4 + idx0
        out(oidx1) = i3i + tr2 + tr3
        out(oidx1 + 1) = i3r + ti2 + ti3
        out(oidx2) = w1r * dr2 - w1i * di2
        out(oidx2 + 1) = w1r * di2 + w1i * dr2
        out(oidx3) = w2r * dr3 - w2i * di3
        out(oidx3 + 1) = w2r * di3 + w2i * dr3
        out(oidx4) = w3r * dr4 - w3i * di4
        out(oidx4 + 1) = w3r * di4 + w3i * dr4
        out(oidx5) = w4r * dr5 - w4i * di5
        out(oidx5 + 1) = w4r * di5 + w4i * dr5

        i += 2
      }
    }
  } /* isign==-1 for forward transform and+1 for backward transform */

  /*----------------------------------------------------------------------
   passfg: Complex FFT's forward/backward processing of general factor;
   isign is +1 for backward and -1 for forward transforms
   ----------------------------------------------------------------------*/
  private[fft] def passfg(nac: Array[Int], ido: Int, ip: Int, l1: Int, idl1: Int, in: Array[Double], in_off: Int, out: Array[Double], out_off: Int, offset: Int, isign: Int): Unit = {
    var idij = 0
    var idlj = 0
    var idot = 0
    var ipph = 0
    var l = 0
    var jc = 0
    var lc = 0
    var idj = 0
    var idl = 0
    var inc = 0
    var idp = 0
    var w1r = .0
    var w1i = .0
    var w2i = .0
    var w2r = .0
    var iw1 = 0
    iw1 = offset
    idot = ido / 2
    ipph = (ip + 1) / 2
    idp = ip * ido
    if (ido >= l1) {
      for (j <- 1 until ipph) {
        jc = ip - j
        val idx1 = j * ido
        val idx2 = jc * ido
        for (k <- 0 until l1) {
          val idx3 = k * ido
          val idx4 = idx3 + idx1 * l1
          val idx5 = idx3 + idx2 * l1
          val idx6 = idx3 * ip
          for (i <- 0 until ido) {
            val oidx1 = out_off + i
            val i1r = in(in_off + i + idx1 + idx6)
            val i2r = in(in_off + i + idx2 + idx6)
            out(oidx1 + idx4) = i1r + i2r
            out(oidx1 + idx5) = i1r - i2r
          }
        }
      }
      for (k <- 0 until l1) {
        val idxt1 = k * ido
        val idxt2 = idxt1 * ip
        for (i <- 0 until ido) {
          out(out_off + i + idxt1) = in(in_off + i + idxt2)
        }
      }
    }
    else {
      for (j <- 1 until ipph) {
        jc = ip - j
        val idxt1 = j * l1 * ido
        val idxt2 = jc * l1 * ido
        val idxt3 = j * ido
        val idxt4 = jc * ido
        for (i <- 0 until ido) {
          for (k <- 0 until l1) {
            val idx1 = k * ido
            val idx2 = idx1 * ip
            val idx3 = out_off + i
            val idx4 = in_off + i
            val i1r = in(idx4 + idxt3 + idx2)
            val i2r = in(idx4 + idxt4 + idx2)
            out(idx3 + idx1 + idxt1) = i1r + i2r
            out(idx3 + idx1 + idxt2) = i1r - i2r
          }
        }
      }
      for (i <- 0 until ido) {
        for (k <- 0 until l1) {
          val idx1 = k * ido
          out(out_off + i + idx1) = in(in_off + i + idx1 * ip)
        }
      }
    }
    idl = 2 - ido
    inc = 0
    val idxt0 = (ip - 1) * idl1
    l = 1
    while ( {
      l < ipph
    }) {
      lc = ip - l
      idl += ido
      val idxt1 = l * idl1
      val idxt2 = lc * idl1
      val idxt3 = idl + iw1
      w1r = wtable(idxt3 - 2)
      w1i = isign * wtable(idxt3 - 1)
      for (ik <- 0 until idl1) {
        val idx1 = in_off + ik
        val idx2 = out_off + ik
        in(idx1 + idxt1) = out(idx2) + w1r * out(idx2 + idl1)
        in(idx1 + idxt2) = w1i * out(idx2 + idxt0)
      }
      idlj = idl
      inc += ido
      for (j <- 2 until ipph) {
        jc = ip - j
        idlj += inc
        if (idlj > idp) idlj -= idp
        val idxt4 = idlj + iw1
        w2r = wtable(idxt4 - 2)
        w2i = isign * wtable(idxt4 - 1)
        val idxt5 = j * idl1
        val idxt6 = jc * idl1
        for (ik <- 0 until idl1) {
          val idx1 = in_off + ik
          val idx2 = out_off + ik
          in(idx1 + idxt1) += w2r * out(idx2 + idxt5)
          in(idx1 + idxt2) += w2i * out(idx2 + idxt6)
        }
      }

      l += 1
    }
    for (j <- 1 until ipph) {
      val idxt1 = j * idl1
      for (ik <- 0 until idl1) {
        val idx1 = out_off + ik
        out(idx1) += out(idx1 + idxt1)
      }
    }
    for (j <- 1 until ipph) {
      jc = ip - j
      val idx1 = j * idl1
      val idx2 = jc * idl1
      var ik = 1
      while ( {
        ik < idl1
      }) {
        val idx3 = out_off + ik
        val idx4 = in_off + ik
        val iidx1 = idx4 + idx1
        val iidx2 = idx4 + idx2
        val i1i = in(iidx1 - 1)
        val i1r = in(iidx1)
        val i2i = in(iidx2 - 1)
        val i2r = in(iidx2)
        val oidx1 = idx3 + idx1
        val oidx2 = idx3 + idx2
        out(oidx1 - 1) = i1i - i2r
        out(oidx2 - 1) = i1i + i2r
        out(oidx1) = i1r + i2i
        out(oidx2) = i1r - i2i

        ik += 2
      }
    }
    nac(0) = 1
    if (ido == 2) return
    nac(0) = 0
    System.arraycopy(out, out_off, in, in_off, idl1)
    val idx0 = l1 * ido
    for (j <- 1 until ip) {
      val idx1 = j * idx0
      for (k <- 0 until l1) {
        val idx2 = k * ido
        val oidx1 = out_off + idx2 + idx1
        val iidx1 = in_off + idx2 + idx1
        in(iidx1) = out(oidx1)
        in(iidx1 + 1) = out(oidx1 + 1)
      }
    }
    if (idot <= l1) {
      idij = 0
      for (j <- 1 until ip) {
        idij += 2
        val idx1 = j * l1 * ido
        var i = 3
        while ( {
          i < ido
        }) {
          idij += 2
          val idx2 = idij + iw1 - 1
          w1r = wtable(idx2 - 1)
          w1i = isign * wtable(idx2)
          val idx3 = in_off + i
          val idx4 = out_off + i
          for (k <- 0 until l1) {
            val idx5 = k * ido + idx1
            val iidx1 = idx3 + idx5
            val oidx1 = idx4 + idx5
            val o1i = out(oidx1 - 1)
            val o1r = out(oidx1)
            in(iidx1 - 1) = w1r * o1i - w1i * o1r
            in(iidx1) = w1r * o1r + w1i * o1i
          }

          i += 2
        }
      }
    }
    else {
      idj = 2 - ido
      for (j <- 1 until ip) {
        idj += ido
        val idx1 = j * l1 * ido
        for (k <- 0 until l1) {
          idij = idj
          val idx3 = k * ido + idx1
          var i = 3
          while ( {
            i < ido
          }) {
            idij += 2
            val idx2 = idij - 1 + iw1
            w1r = wtable(idx2 - 1)
            w1i = isign * wtable(idx2)
            val iidx1 = in_off + i + idx3
            val oidx1 = out_off + i + idx3
            val o1i = out(oidx1 - 1)
            val o1r = out(oidx1)
            in(iidx1 - 1) = w1r * o1i - w1i * o1r
            in(iidx1) = w1r * o1r + w1i * o1i

            i += 2
          }
        }
      }
    }
  }
}
