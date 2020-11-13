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

import java.lang.Math.min

import de.sciss.transform4s.utils.ConcurrencyUtils.executionContext
import de.sciss.transform4s.utils.{CommonUtils, ConcurrencyUtils}

import scala.concurrent.Future

/**
 * Computes 2D Discrete Fourier Transform (DFT) of complex and real, double
 * precision data. The sizes of both dimensions can be arbitrary numbers. This
 * is a parallel implementation of split-radix and mixed-radix algorithms
 * optimized for SMP systems. <br>
 * <br>
 * Part of the code is derived from General Purpose FFT Package written by Takuya Ooura
 * (http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html)
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
object DoubleFFT_2D {
  def apply(rows: Int, columns: Int): DoubleFFT_2D = {
    if (rows <= 1 || columns <= 1) throw new IllegalArgumentException("rows and columns must be greater than 1")

    val useThreads    = rows * columns >= CommonUtils.threadsBeginN_2D
    val isPowerOfTwo  = CommonUtils.isPowerOf2(rows) && CommonUtils.isPowerOf2(columns)
    val fftRows       = DoubleFFT_1D(rows)
    val fftColumns    = if (rows == columns) fftRows else DoubleFFT_1D(columns)

    new DoubleFFT_2D(
      rows        = rows        ,
      columns0    = columns     ,
      fftColumns  = fftColumns  ,
      fftRows     = fftRows     ,
      isPowerOfTwo= isPowerOfTwo,
      useThreads  = useThreads  ,
    )
  }
}
class DoubleFFT_2D private(
                            rows        : Int         ,
                            columns0    : Int         ,
                            fftColumns  : DoubleFFT_1D,
                            fftRows     : DoubleFFT_1D,
                            isPowerOfTwo: Boolean     ,
                            useThreads  : Boolean     ,
                          ) {

  private var columns: Int = columns0

  /**
   * Computes 2D forward DFT of complex data leaving the result in
   * <code>a</code>. The data is stored in 1D array in row-major order.
   * Complex number is stored as two double values in sequence: the real and
   * imaginary part, i.e. the input array must be of size rows*2*columns. The
   * physical layout of the input data has to be as follows:<br>
   *
   * <pre>
   * a[k1*2*columns+2*k2] = Re[k1][k2],
   * a[k1*2*columns+2*k2+1] = Im[k1][k2], 0&lt;=k1&lt;rows, 0&lt;=k2&lt;columns,
   * </pre>
   *
   * @param a
   * data to transform
   */
  def complexForward(a: Array[Double]): Unit = {
    val nthreads: Int = ConcurrencyUtils.numThreads
    if (isPowerOfTwo) {
      columns = 2 * columns
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(0, -1, a, scale = true)
        cdft2d_subth(-1, a, scale = true)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexForward(a, r * columns)
        }
        cdft2d_sub(-1, a, scale = true)
      }
      columns = columns / 2
    }
    else {
      val rowStride: Int = 2 * columns
      if ((nthreads > 1) && useThreads && (rows >= nthreads) && (columns >= nthreads)) {
        val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
        var p: Int = rows / nthreads
        for (l <- 0 until nthreads) {
          val firstRow: Int = l * p
          val lastRow: Int = if (l == (nthreads - 1)) {
            rows
          }
          else {
            firstRow + p
          }
          futures(l) = Future {
            for (r <- firstRow until lastRow) {
              fftColumns.complexForward(a, r * rowStride)
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
        p = columns / nthreads
        for (l <- 0 until nthreads) {
          val firstColumn: Int = l * p
          val lastColumn: Int = if (l == (nthreads - 1)) {
            columns
          }
          else {
            firstColumn + p
          }
          futures(l) = Future {
            val temp: Array[Double] = new Array[Double](2 * rows)
            for (c <- firstColumn until lastColumn) {
              val idx0: Int = 2 * c
              for (r <- 0 until rows) {
                val idx1: Int = 2 * r
                val idx2: Int = r * rowStride + idx0
                temp(idx1) = a(idx2)
                temp(idx1 + 1) = a(idx2 + 1)
              }
              fftRows.complexForward(temp)
              for (r <- 0 until rows) {
                val idx1: Int = 2 * r
                val idx2: Int = r * rowStride + idx0
                a(idx2) = temp(idx1)
                a(idx2 + 1) = temp(idx1 + 1)
              }
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexForward(a, r * rowStride)
        }
        val temp: Array[Double] = new Array[Double](2 * rows)
        for (c <- 0 until columns) {
          val idx0: Int = 2 * c
          for (r <- 0 until rows) {
            val idx1: Int = 2 * r
            val idx2: Int = r * rowStride + idx0
            temp(idx1) = a(idx2)
            temp(idx1 + 1) = a(idx2 + 1)
          }
          fftRows.complexForward(temp)
          for (r <- 0 until rows) {
            val idx1: Int = 2 * r
            val idx2: Int = r * rowStride + idx0
            a(idx2) = temp(idx1)
            a(idx2 + 1) = temp(idx1 + 1)
          }
        }
      }
    }
  }

  /**
   * Computes 2D forward DFT of complex data leaving the result in
   * <code>a</code>. The data is stored in 2D array. Complex data is
   * represented by 2 double values in sequence: the real and imaginary part,
   * i.e. the input array must be of size rows by 2*columns. The physical
   * layout of the input data has to be as follows:<br>
   *
   * <pre>
   * a[k1][2*k2] = Re[k1][k2],
   * a[k1][2*k2+1] = Im[k1][k2], 0&lt;=k1&lt;rows, 0&lt;=k2&lt;columns,
   * </pre>
   *
   * @param a
   * data to transform
   */
  def complexForward(a: Array[Array[Double]]): Unit = {
    val nthreads: Int = ConcurrencyUtils.numThreads
    if (isPowerOfTwo) {
      columns = 2 * columns
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(0, -1, a, scale = true)
        cdft2d_subth(-1, a, scale = true)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexForward(a(r))
        }
        cdft2d_sub(-1, a, scale = true)
      }
      columns = columns / 2
    }
    else {
      if ((nthreads > 1) && useThreads && (rows >= nthreads) && (columns >= nthreads)) {
        val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
        var p: Int = rows / nthreads
        for (l <- 0 until nthreads) {
          val firstRow: Int = l * p
          val lastRow: Int = if (l == (nthreads - 1)) {
            rows
          }
          else {
            firstRow + p
          }
          futures(l) = Future {
            for (r <- firstRow until lastRow) {
              fftColumns.complexForward(a(r))
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
        p = columns / nthreads
        for (l <- 0 until nthreads) {
          val firstColumn: Int = l * p
          val lastColumn: Int = if (l == (nthreads - 1)) {
            columns
          }
          else {
            firstColumn + p
          }
          futures(l) = Future {
            val temp: Array[Double] = new Array[Double](2 * rows)
            for (c <- firstColumn until lastColumn) {
              val idx1: Int = 2 * c
              for (r <- 0 until rows) {
                val idx2: Int = 2 * r
                temp(idx2) = a(r)(idx1)
                temp(idx2 + 1) = a(r)(idx1 + 1)
              }
              fftRows.complexForward(temp)
              for (r <- 0 until rows) {
                val idx2: Int = 2 * r
                a(r)(idx1) = temp(idx2)
                a(r)(idx1 + 1) = temp(idx2 + 1)
              }
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexForward(a(r))
        }
        val temp: Array[Double] = new Array[Double](2 * rows)
        for (c <- 0 until columns) {
          val idx1: Int = 2 * c
          for (r <- 0 until rows) {
            val idx2: Int = 2 * r
            temp(idx2) = a(r)(idx1)
            temp(idx2 + 1) = a(r)(idx1 + 1)
          }
          fftRows.complexForward(temp)
          for (r <- 0 until rows) {
            val idx2: Int = 2 * r
            a(r)(idx1) = temp(idx2)
            a(r)(idx1 + 1) = temp(idx2 + 1)
          }
        }
      }
    }
  }

  /**
   * Computes 2D inverse DFT of complex data leaving the result in
   * <code>a</code>. The data is stored in 1D array in row-major order.
   * Complex number is stored as two double values in sequence: the real and
   * imaginary part, i.e. the input array must be of size rows*2*columns. The
   * physical layout of the input data has to be as follows:<br>
   *
   * <pre>
   * a[k1*2*columns+2*k2] = Re[k1][k2],
   * a[k1*2*columns+2*k2+1] = Im[k1][k2], 0&lt;=k1&lt;rows, 0&lt;=k2&lt;columns,
   * </pre>
   *
   * @param a
   * data to transform
   * @param scale
   * if true then scaling is performed
   *
   */
  def complexInverse(a: Array[Double], scale: Boolean): Unit = {
    val nthreads: Int = ConcurrencyUtils.numThreads
    if (isPowerOfTwo) {
      columns = 2 * columns
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(0, 1, a, scale)
        cdft2d_subth(1, a, scale)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexInverse(a, r * columns, scale)
        }
        cdft2d_sub(1, a, scale)
      }
      columns = columns / 2
    }
    else {
      val rowspan: Int = 2 * columns
      if ((nthreads > 1) && useThreads && (rows >= nthreads) && (columns >= nthreads)) {
        val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
        var p: Int = rows / nthreads
        for (l <- 0 until nthreads) {
          val firstRow: Int = l * p
          val lastRow: Int = if (l == (nthreads - 1)) {
            rows
          }
          else {
            firstRow + p
          }
          futures(l) = Future {
            for (r <- firstRow until lastRow) {
              fftColumns.complexInverse(a, r * rowspan, scale)
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
        p = columns / nthreads
        for (l <- 0 until nthreads) {
          val firstColumn: Int = l * p
          val lastColumn: Int = if (l == (nthreads - 1)) {
            columns
          }
          else {
            firstColumn + p
          }
          futures(l) = Future {
            val temp: Array[Double] = new Array[Double](2 * rows)
            for (c <- firstColumn until lastColumn) {
              val idx1: Int = 2 * c
              for (r <- 0 until rows) {
                val idx2: Int = 2 * r
                val idx3: Int = r * rowspan + idx1
                temp(idx2) = a(idx3)
                temp(idx2 + 1) = a(idx3 + 1)
              }
              fftRows.complexInverse(temp, scale)
              for (r <- 0 until rows) {
                val idx2: Int = 2 * r
                val idx3: Int = r * rowspan + idx1
                a(idx3) = temp(idx2)
                a(idx3 + 1) = temp(idx2 + 1)
              }
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexInverse(a, r * rowspan, scale)
        }
        val temp: Array[Double] = new Array[Double](2 * rows)
        for (c <- 0 until columns) {
          val idx1: Int = 2 * c
          for (r <- 0 until rows) {
            val idx2: Int = 2 * r
            val idx3: Int = r * rowspan + idx1
            temp(idx2) = a(idx3)
            temp(idx2 + 1) = a(idx3 + 1)
          }
          fftRows.complexInverse(temp, scale)
          for (r <- 0 until rows) {
            val idx2: Int = 2 * r
            val idx3: Int = r * rowspan + idx1
            a(idx3) = temp(idx2)
            a(idx3 + 1) = temp(idx2 + 1)
          }
        }
      }
    }
  }

  /**
   * Computes 2D inverse DFT of complex data leaving the result in
   * <code>a</code>. The data is stored in 2D array. Complex data is
   * represented by 2 double values in sequence: the real and imaginary part,
   * i.e. the input array must be of size rows by 2*columns. The physical
   * layout of the input data has to be as follows:<br>
   *
   * <pre>
   * a[k1][2*k2] = Re[k1][k2],
   * a[k1][2*k2+1] = Im[k1][k2], 0&lt;=k1&lt;rows, 0&lt;=k2&lt;columns,
   * </pre>
   *
   * @param a
   * data to transform
   * @param scale
   * if true then scaling is performed
   *
   */
  def complexInverse(a: Array[Array[Double]], scale: Boolean): Unit = {
    val nthreads: Int = ConcurrencyUtils.numThreads
    if (isPowerOfTwo) {
      columns = 2 * columns
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(0, 1, a, scale)
        cdft2d_subth(1, a, scale)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexInverse(a(r), scale)
        }
        cdft2d_sub(1, a, scale)
      }
      columns = columns / 2
    }
    else {
      if ((nthreads > 1) && useThreads && (rows >= nthreads) && (columns >= nthreads)) {
        val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
        var p: Int = rows / nthreads
        for (l <- 0 until nthreads) {
          val firstRow: Int = l * p
          val lastRow: Int = if (l == (nthreads - 1)) {
            rows
          }
          else {
            firstRow + p
          }
          futures(l) = Future {
            for (r <- firstRow until lastRow) {
              fftColumns.complexInverse(a(r), scale)
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
        p = columns / nthreads
        for (l <- 0 until nthreads) {
          val firstColumn: Int = l * p
          val lastColumn: Int = if (l == (nthreads - 1)) {
            columns
          }
          else {
            firstColumn + p
          }
          futures(l) = Future {
            val temp: Array[Double] = new Array[Double](2 * rows)
            for (c <- firstColumn until lastColumn) {
              val idx1: Int = 2 * c
              for (r <- 0 until rows) {
                val idx2: Int = 2 * r
                temp(idx2) = a(r)(idx1)
                temp(idx2 + 1) = a(r)(idx1 + 1)
              }
              fftRows.complexInverse(temp, scale)
              for (r <- 0 until rows) {
                val idx2: Int = 2 * r
                a(r)(idx1) = temp(idx2)
                a(r)(idx1 + 1) = temp(idx2 + 1)
              }
            }
          }
        }
        ConcurrencyUtils.waitForCompletion(futures)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.complexInverse(a(r), scale)
        }
        val temp: Array[Double] = new Array[Double](2 * rows)
        for (c <- 0 until columns) {
          val idx1: Int = 2 * c
          for (r <- 0 until rows) {
            val idx2: Int = 2 * r
            temp(idx2) = a(r)(idx1)
            temp(idx2 + 1) = a(r)(idx1 + 1)
          }
          fftRows.complexInverse(temp, scale)
          for (r <- 0 until rows) {
            val idx2: Int = 2 * r
            a(r)(idx1) = temp(idx2)
            a(r)(idx1 + 1) = temp(idx2 + 1)
          }
        }
      }
    }
  }

  /**
   * Computes 2D forward DFT of real data leaving the result in <code>a</code>
   * . This method only works when the sizes of both dimensions are
   * power-of-two numbers. The physical layout of the output data is as
   * follows:
   *
   * <pre>
   * a[k1*columns+2*k2] = Re[k1][k2] = Re[rows-k1][columns-k2],
   * a[k1*columns+2*k2+1] = Im[k1][k2] = -Im[rows-k1][columns-k2],
   * 0&lt;k1&lt;rows, 0&lt;k2&lt;columns/2,
   * a[2*k2] = Re[0][k2] = Re[0][columns-k2],
   * a[2*k2+1] = Im[0][k2] = -Im[0][columns-k2],
   * 0&lt;k2&lt;columns/2,
   * a[k1*columns] = Re[k1][0] = Re[rows-k1][0],
   * a[k1*columns+1] = Im[k1][0] = -Im[rows-k1][0],
   * a[(rows-k1)*columns+1] = Re[k1][columns/2] = Re[rows-k1][columns/2],
   * a[(rows-k1)*columns] = -Im[k1][columns/2] = Im[rows-k1][columns/2],
   * 0&lt;k1&lt;rows/2,
   * a[0] = Re[0][0],
   * a[1] = Re[0][columns/2],
   * a[(rows/2)*columns] = Re[rows/2][0],
   * a[(rows/2)*columns+1] = Re[rows/2][columns/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * forward transform, use <code>realForwardFull</code>. To get back the
   * original data, use <code>realInverse</code> on the output of this method.
   *
   * @param a
   * data to transform
   */
  def realForward(a: Array[Double]): Unit = {
    if (!isPowerOfTwo) {
      throw new IllegalArgumentException("rows and columns must be power of two numbers")
    }
    else {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(1, 1, a, scale = true)
        cdft2d_subth(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.realForward(a, r * columns)
        }
        cdft2d_sub(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
    }
  }

  /**
   * Computes 2D forward DFT of real data leaving the result in <code>a</code>
   * . This method only works when the sizes of both dimensions are
   * power-of-two numbers. The physical layout of the output data is as
   * follows:
   *
   * <pre>
   * a[k1][2*k2] = Re[k1][k2] = Re[rows-k1][columns-k2],
   * a[k1][2*k2+1] = Im[k1][k2] = -Im[rows-k1][columns-k2],
   * 0&lt;k1&lt;rows, 0&lt;k2&lt;columns/2,
   * a[0][2*k2] = Re[0][k2] = Re[0][columns-k2],
   * a[0][2*k2+1] = Im[0][k2] = -Im[0][columns-k2],
   * 0&lt;k2&lt;columns/2,
   * a[k1][0] = Re[k1][0] = Re[rows-k1][0],
   * a[k1][1] = Im[k1][0] = -Im[rows-k1][0],
   * a[rows-k1][1] = Re[k1][columns/2] = Re[rows-k1][columns/2],
   * a[rows-k1][0] = -Im[k1][columns/2] = Im[rows-k1][columns/2],
   * 0&lt;k1&lt;rows/2,
   * a[0][0] = Re[0][0],
   * a[0][1] = Re[0][columns/2],
   * a[rows/2][0] = Re[rows/2][0],
   * a[rows/2][1] = Re[rows/2][columns/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * forward transform, use <code>realForwardFull</code>. To get back the
   * original data, use <code>realInverse</code> on the output of this method.
   *
   * @param a
   * data to transform
   */
  def realForward(a: Array[Array[Double]]): Unit = {
    if (!isPowerOfTwo) {
      throw new IllegalArgumentException("rows and columns must be power of two numbers")
    }
    else {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(1, 1, a, scale = true)
        cdft2d_subth(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.realForward(a(r))
        }
        cdft2d_sub(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
    }
  }

  /**
   * Computes 2D forward DFT of real data leaving the result in <code>a</code>
   * . This method computes full real forward transform, i.e. you will get the
   * same result as from <code>complexForward</code> called with all imaginary
   * part equal 0. Because the result is stored in <code>a</code>, the input
   * array must be of size rows*2*columns, with only the first rows*columns
   * elements filled with real data. To get back the original data, use
   * <code>complexInverse</code> on the output of this method.
   *
   * @param a
   * data to transform
   */
  def realForwardFull(a: Array[Double]): Unit = {
    if (isPowerOfTwo) {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(1, 1, a, scale = true)
        cdft2d_subth(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.realForward(a, r * columns)
        }
        cdft2d_sub(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
      fillSymmetric(a)
    }
    else {
      mixedRadixRealForwardFull(a)
    }
  }

  /**
   * Computes 2D forward DFT of real data leaving the result in <code>a</code>
   * . This method computes full real forward transform, i.e. you will get the
   * same result as from <code>complexForward</code> called with all imaginary
   * part equal 0. Because the result is stored in <code>a</code>, the input
   * array must be of size rows by 2*columns, with only the first rows by
   * columns elements filled with real data. To get back the original data,
   * use <code>complexInverse</code> on the output of this method.
   *
   * @param a
   * data to transform
   */
  def realForwardFull(a: Array[Array[Double]]): Unit = {
    if (isPowerOfTwo) {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth1(1, 1, a, scale = true)
        cdft2d_subth(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.realForward(a(r))
        }
        cdft2d_sub(-1, a, scale = true)
        rdft2d_sub(1, a)
      }
      fillSymmetric(a)
    }
    else {
      mixedRadixRealForwardFull(a)
    }
  }

  /**
   * Computes 2D inverse DFT of real data leaving the result in <code>a</code>
   * . This method only works when the sizes of both dimensions are
   * power-of-two numbers. The physical layout of the input data has to be as
   * follows:
   *
   * <pre>
   * a[k1*columns+2*k2] = Re[k1][k2] = Re[rows-k1][columns-k2],
   * a[k1*columns+2*k2+1] = Im[k1][k2] = -Im[rows-k1][columns-k2],
   * 0&lt;k1&lt;rows, 0&lt;k2&lt;columns/2,
   * a[2*k2] = Re[0][k2] = Re[0][columns-k2],
   * a[2*k2+1] = Im[0][k2] = -Im[0][columns-k2],
   * 0&lt;k2&lt;columns/2,
   * a[k1*columns] = Re[k1][0] = Re[rows-k1][0],
   * a[k1*columns+1] = Im[k1][0] = -Im[rows-k1][0],
   * a[(rows-k1)*columns+1] = Re[k1][columns/2] = Re[rows-k1][columns/2],
   * a[(rows-k1)*columns] = -Im[k1][columns/2] = Im[rows-k1][columns/2],
   * 0&lt;k1&lt;rows/2,
   * a[0] = Re[0][0],
   * a[1] = Re[0][columns/2],
   * a[(rows/2)*columns] = Re[rows/2][0],
   * a[(rows/2)*columns+1] = Re[rows/2][columns/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * inverse transform, use <code>realInverseFull</code>.
   *
   * @param a
   * data to transform
   * @param scale
   * if true then scaling is performed
   */
  def realInverse(a: Array[Double], scale: Boolean): Unit = {
    if (!isPowerOfTwo) {
      throw new IllegalArgumentException("rows and columns must be power of two numbers")
    }
    else {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        rdft2d_sub(-1, a)
        cdft2d_subth(1, a, scale)
        xdft2d0_subth1(1, -1, a, scale)
      }
      else {
        rdft2d_sub(-1, a)
        cdft2d_sub(1, a, scale)
        for (r <- 0 until rows) {
          fftColumns.realInverse(a, r * columns, scale)
        }
      }
    }
  }

  /**
   * Computes 2D inverse DFT of real data leaving the result in <code>a</code>
   * . This method only works when the sizes of both dimensions are
   * power-of-two numbers. The physical layout of the input data has to be as
   * follows:
   *
   * <pre>
   * a[k1][2*k2] = Re[k1][k2] = Re[rows-k1][columns-k2],
   * a[k1][2*k2+1] = Im[k1][k2] = -Im[rows-k1][columns-k2],
   * 0&lt;k1&lt;rows, 0&lt;k2&lt;columns/2,
   * a[0][2*k2] = Re[0][k2] = Re[0][columns-k2],
   * a[0][2*k2+1] = Im[0][k2] = -Im[0][columns-k2],
   * 0&lt;k2&lt;columns/2,
   * a[k1][0] = Re[k1][0] = Re[rows-k1][0],
   * a[k1][1] = Im[k1][0] = -Im[rows-k1][0],
   * a[rows-k1][1] = Re[k1][columns/2] = Re[rows-k1][columns/2],
   * a[rows-k1][0] = -Im[k1][columns/2] = Im[rows-k1][columns/2],
   * 0&lt;k1&lt;rows/2,
   * a[0][0] = Re[0][0],
   * a[0][1] = Re[0][columns/2],
   * a[rows/2][0] = Re[rows/2][0],
   * a[rows/2][1] = Re[rows/2][columns/2]
   * </pre>
   *
   * This method computes only half of the elements of the real transform. The
   * other half satisfies the symmetry condition. If you want the full real
   * inverse transform, use <code>realInverseFull</code>.
   *
   * @param a
   * data to transform
   * @param scale
   * if true then scaling is performed
   */
  def realInverse(a: Array[Array[Double]], scale: Boolean): Unit = {
    if (!isPowerOfTwo) {
      throw new IllegalArgumentException("rows and columns must be power of two numbers")
    }
    else {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        rdft2d_sub(-1, a)
        cdft2d_subth(1, a, scale)
        xdft2d0_subth1(1, -1, a, scale)
      }
      else {
        rdft2d_sub(-1, a)
        cdft2d_sub(1, a, scale)
        for (r <- 0 until rows) {
          fftColumns.realInverse(a(r), scale)
        }
      }
    }
  }

  /**
   * Computes 2D inverse DFT of real data leaving the result in <code>a</code>
   * . This method computes full real inverse transform, i.e. you will get the
   * same result as from <code>complexInverse</code> called with all imaginary
   * part equal 0. Because the result is stored in <code>a</code>, the input
   * array must be of size rows*2*columns, with only the first rows*columns
   * elements filled with real data.
   *
   * @param a
   * data to transform
   * @param scale
   * if true then scaling is performed
   */
  def realInverseFull(a: Array[Double], scale: Boolean): Unit = {
    if (isPowerOfTwo) {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth2(1, -1, a, scale)
        cdft2d_subth(1, a, scale)
        rdft2d_sub(1, a)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.realInverse2(a, r * columns, scale)
        }
        cdft2d_sub(1, a, scale)
        rdft2d_sub(1, a)
      }
      fillSymmetric(a)
    }
    else {
      mixedRadixRealInverseFull(a, scale)
    }
  }

  /**
   * Computes 2D inverse DFT of real data leaving the result in <code>a</code>
   * . This method computes full real inverse transform, i.e. you will get the
   * same result as from <code>complexInverse</code> called with all imaginary
   * part equal 0. Because the result is stored in <code>a</code>, the input
   * array must be of size rows by 2*columns, with only the first rows by
   * columns elements filled with real data.
   *
   * @param a
   * data to transform
   * @param scale
   * if true then scaling is performed
   */
  def realInverseFull(a: Array[Array[Double]], scale: Boolean): Unit = {
    if (isPowerOfTwo) {
      val nthreads: Int = ConcurrencyUtils.numThreads
      if ((nthreads > 1) && useThreads) {
        xdft2d0_subth2(1, -1, a, scale)
        cdft2d_subth(1, a, scale)
        rdft2d_sub(1, a)
      }
      else {
        for (r <- 0 until rows) {
          fftColumns.realInverse2(a(r), 0, scale)
        }
        cdft2d_sub(1, a, scale)
        rdft2d_sub(1, a)
      }
      fillSymmetric(a)
    }
    else {
      mixedRadixRealInverseFull(a, scale)
    }
  }

  private def mixedRadixRealForwardFull(a: Array[Array[Double]]): Unit = {
    val n2d2: Int = columns / 2 + 1
    val temp: Array[Array[Double]] = Array.ofDim[Double](n2d2, 2 * rows)
    val nthreads: Int = ConcurrencyUtils.numThreads
    if ((nthreads > 1) && useThreads && (rows >= nthreads) && (n2d2 - 2 >= nthreads)) {
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var p: Int = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future  {
          for (i <- firstRow until lastRow) {
            fftColumns.realForward(a(i))
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (r <- 0 until rows) {
        temp(0)(r) = a(r)(0) //first column is always real

      }
      fftRows.realForwardFull(temp(0))
      p = (n2d2 - 2) / nthreads
      for (l <- 0 until nthreads) {
        val firstColumn: Int = 1 + l * p
        val lastColumn: Int = if (l == (nthreads - 1)) {
          n2d2 - 1
        }
        else {
          firstColumn + p
        }
        futures(l) = Future {
          for (c <- firstColumn until lastColumn) {
            val idx2: Int = 2 * c
            for (r <- 0 until rows) {
              val idx1: Int = 2 * r
              temp(c)(idx1) = a(r)(idx2)
              temp(c)(idx1 + 1) = a(r)(idx2 + 1)
            }
            fftRows.complexForward(temp(c))
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r)(1)
          //imaginary part = 0;
        }
        fftRows.realForwardFull(temp(n2d2 - 1))
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = n2d2 - 1
          temp(idx2)(idx1) = a(r)(2 * idx2)
          temp(idx2)(idx1 + 1) = a(r)(1)
        }
        fftRows.complexForward(temp(n2d2 - 1))
      }
      p = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx1: Int = 2 * r
            for (c <- 0 until n2d2) {
              val idx2: Int = 2 * c
              a(r)(idx2) = temp(c)(idx1)
              a(r)(idx2 + 1) = temp(c)(idx1 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (l <- 0 until nthreads) {
        val firstRow: Int = 1 + l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx3: Int = rows - r
            for (c <- n2d2 until columns) {
              val idx1: Int = 2 * c
              val idx2: Int = 2 * (columns - c)
              a(0)(idx1) = a(0)(idx2)
              a(0)(idx1 + 1) = -a(0)(idx2 + 1)
              a(r)(idx1) = a(idx3)(idx2)
              a(r)(idx1 + 1) = -a(idx3)(idx2 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (r <- 0 until rows) {
        fftColumns.realForward(a(r))
      }
      for (r <- 0 until rows) {
        temp(0)(r) = a(r)(0)
      }
      fftRows.realForwardFull(temp(0))
      for (c <- 1 until n2d2 - 1) {
        val idx2: Int = 2 * c
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          temp(c)(idx1) = a(r)(idx2)
          temp(c)(idx1 + 1) = a(r)(idx2 + 1)
        }
        fftRows.complexForward(temp(c))
      }
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r)(1)
        }
        fftRows.realForwardFull(temp(n2d2 - 1))
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = n2d2 - 1
          temp(idx2)(idx1) = a(r)(2 * idx2)
          temp(idx2)(idx1 + 1) = a(r)(1)
        }
        fftRows.complexForward(temp(n2d2 - 1))
      }
      for (r <- 0 until rows) {
        val idx1: Int = 2 * r
        for (c <- 0 until n2d2) {
          val idx2: Int = 2 * c
          a(r)(idx2) = temp(c)(idx1)
          a(r)(idx2 + 1) = temp(c)(idx1 + 1)
        }
      }
      //fill symmetric
      for (r <- 1 until rows) {
        val idx3: Int = rows - r
        for (c <- n2d2 until columns) {
          val idx1: Int = 2 * c
          val idx2: Int = 2 * (columns - c)
          a(0)(idx1) = a(0)(idx2)
          a(0)(idx1 + 1) = -a(0)(idx2 + 1)
          a(r)(idx1) = a(idx3)(idx2)
          a(r)(idx1 + 1) = -a(idx3)(idx2 + 1)
        }
      }
    }
  }

  private def mixedRadixRealForwardFull(a: Array[Double]): Unit = {
    val rowStride: Int = 2 * columns
    val n2d2: Int = columns / 2 + 1
    val temp: Array[Array[Double]] = Array.ofDim[Double](n2d2, 2 * rows)
    val nthreads: Int = ConcurrencyUtils.numThreads
    if ((nthreads > 1) && useThreads && (rows >= nthreads) && (n2d2 - 2 >= nthreads)) {
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var p: Int = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (i <- firstRow until lastRow) {
            fftColumns.realForward(a, i * columns)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (r <- 0 until rows) {
        temp(0)(r) = a(r * columns)
      }
      fftRows.realForwardFull(temp(0))
      p = (n2d2 - 2) / nthreads
      for (l <- 0 until nthreads) {
        val firstColumn: Int = 1 + l * p
        val lastColumn: Int = if (l == (nthreads - 1)) {
          n2d2 - 1
        }
        else {
          firstColumn + p
        }
        futures(l) = Future {
          for (c <- firstColumn until lastColumn) {
            val idx0: Int = 2 * c
            for (r <- 0 until rows) {
              val idx1: Int = 2 * r
              val idx2: Int = r * columns + idx0
              temp(c)(idx1) = a(idx2)
              temp(c)(idx1 + 1) = a(idx2 + 1)
            }
            fftRows.complexForward(temp(c))
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r * columns + 1)
        }
        fftRows.realForwardFull(temp(n2d2 - 1))
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = r * columns
          val idx3: Int = n2d2 - 1
          temp(idx3)(idx1) = a(idx2 + 2 * idx3)
          temp(idx3)(idx1 + 1) = a(idx2 + 1)
        }
        fftRows.complexForward(temp(n2d2 - 1))
      }
      p = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx1: Int = 2 * r
            for (c <- 0 until n2d2) {
              val idx0: Int = 2 * c
              val idx2: Int = r * rowStride + idx0
              a(idx2) = temp(c)(idx1)
              a(idx2 + 1) = temp(c)(idx1 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (l <- 0 until nthreads) {
        val firstRow: Int = 1 + l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx5: Int = r * rowStride
            val idx6: Int = (rows - r + 1) * rowStride
            for (c <- n2d2 until columns) {
              val idx1: Int = 2 * c
              val idx2: Int = 2 * (columns - c)
              a(idx1) = a(idx2)
              a(idx1 + 1) = -a(idx2 + 1)
              val idx3: Int = idx5 + idx1
              val idx4: Int = idx6 - idx1
              a(idx3) = a(idx4)
              a(idx3 + 1) = -a(idx4 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (r <- 0 until rows) {
        fftColumns.realForward(a, r * columns)
      }
      for (r <- 0 until rows) {
        temp(0)(r) = a(r * columns)
      }
      fftRows.realForwardFull(temp(0))
      for (c <- 1 until n2d2 - 1) {
        val idx0: Int = 2 * c
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = r * columns + idx0
          temp(c)(idx1) = a(idx2)
          temp(c)(idx1 + 1) = a(idx2 + 1)
        }
        fftRows.complexForward(temp(c))
      }
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r * columns + 1)
        }
        fftRows.realForwardFull(temp(n2d2 - 1))
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = r * columns
          val idx3: Int = n2d2 - 1
          temp(idx3)(idx1) = a(idx2 + 2 * idx3)
          temp(idx3)(idx1 + 1) = a(idx2 + 1)
        }
        fftRows.complexForward(temp(n2d2 - 1))
      }
      for (r <- 0 until rows) {
        val idx1: Int = 2 * r
        for (c <- 0 until n2d2) {
          val idx0: Int = 2 * c
          val idx2: Int = r * rowStride + idx0
          a(idx2) = temp(c)(idx1)
          a(idx2 + 1) = temp(c)(idx1 + 1)
        }
      }
      for (r <- 1 until rows) {
        val idx5: Int = r * rowStride
        val idx6: Int = (rows - r + 1) * rowStride
        for (c <- n2d2 until columns) {
          val idx1: Int = 2 * c
          val idx2: Int = 2 * (columns - c)
          a(idx1) = a(idx2)
          a(idx1 + 1) = -a(idx2 + 1)
          val idx3: Int = idx5 + idx1
          val idx4: Int = idx6 - idx1
          a(idx3) = a(idx4)
          a(idx3 + 1) = -a(idx4 + 1)
        }
      }
    }
  }

  private def mixedRadixRealInverseFull(a: Array[Array[Double]], scale: Boolean): Unit = {
    val n2d2: Int = columns / 2 + 1
    val temp: Array[Array[Double]] = Array.ofDim[Double](n2d2, 2 * rows)
    val nthreads: Int = ConcurrencyUtils.numThreads
    if ((nthreads > 1) && useThreads && (rows >= nthreads) && (n2d2 - 2 >= nthreads)) {
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var p: Int = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (i <- firstRow until lastRow) {
            fftColumns.realInverse2(a(i), 0, scale)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (r <- 0 until rows) {
        temp(0)(r) = a(r)(0)
      }
      fftRows.realInverseFull(temp(0), scale)
      p = (n2d2 - 2) / nthreads
      for (l <- 0 until nthreads) {
        val firstColumn: Int = 1 + l * p
        val lastColumn: Int = if (l == (nthreads - 1)) {
          n2d2 - 1
        }
        else {
          firstColumn + p
        }
        futures(l) = Future {
          for (c <- firstColumn until lastColumn) {
            val idx2: Int = 2 * c
            for (r <- 0 until rows) {
              val idx1: Int = 2 * r
              temp(c)(idx1) = a(r)(idx2)
              temp(c)(idx1 + 1) = a(r)(idx2 + 1)
            }
            fftRows.complexInverse(temp(c), scale)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r)(1)
        }
        fftRows.realInverseFull(temp(n2d2 - 1), scale)
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = n2d2 - 1
          temp(idx2)(idx1) = a(r)(2 * idx2)
          temp(idx2)(idx1 + 1) = a(r)(1)
        }
        fftRows.complexInverse(temp(n2d2 - 1), scale)
      }
      p = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx1: Int = 2 * r
            for (c <- 0 until n2d2) {
              val idx2: Int = 2 * c
              a(r)(idx2) = temp(c)(idx1)
              a(r)(idx2 + 1) = temp(c)(idx1 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (l <- 0 until nthreads) {
        val firstRow: Int = 1 + l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx3: Int = rows - r
            for (c <- n2d2 until columns) {
              val idx1: Int = 2 * c
              val idx2: Int = 2 * (columns - c)
              a(0)(idx1) = a(0)(idx2)
              a(0)(idx1 + 1) = -a(0)(idx2 + 1)
              a(r)(idx1) = a(idx3)(idx2)
              a(r)(idx1 + 1) = -a(idx3)(idx2 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (r <- 0 until rows) {
        fftColumns.realInverse2(a(r), 0, scale)
      }
      for (r <- 0 until rows) {
        temp(0)(r) = a(r)(0)
      }
      fftRows.realInverseFull(temp(0), scale)
      for (c <- 1 until n2d2 - 1) {
        val idx2: Int = 2 * c
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          temp(c)(idx1) = a(r)(idx2)
          temp(c)(idx1 + 1) = a(r)(idx2 + 1)
        }
        fftRows.complexInverse(temp(c), scale)
      }
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r)(1)
        }
        fftRows.realInverseFull(temp(n2d2 - 1), scale)
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = n2d2 - 1
          temp(idx2)(idx1) = a(r)(2 * idx2)
          temp(idx2)(idx1 + 1) = a(r)(1)
        }
        fftRows.complexInverse(temp(n2d2 - 1), scale)
      }
      for (r <- 0 until rows) {
        val idx1: Int = 2 * r
        for (c <- 0 until n2d2) {
          val idx2: Int = 2 * c
          a(r)(idx2) = temp(c)(idx1)
          a(r)(idx2 + 1) = temp(c)(idx1 + 1)
        }
      }
      for (r <- 1 until rows) {
        val idx3: Int = rows - r
        for (c <- n2d2 until columns) {
          val idx1: Int = 2 * c
          val idx2: Int = 2 * (columns - c)
          a(0)(idx1) = a(0)(idx2)
          a(0)(idx1 + 1) = -a(0)(idx2 + 1)
          a(r)(idx1) = a(idx3)(idx2)
          a(r)(idx1 + 1) = -a(idx3)(idx2 + 1)
        }
      }
    }
  }

  private def mixedRadixRealInverseFull(a: Array[Double], scale: Boolean): Unit = {
    val rowStride: Int = 2 * columns
    val n2d2: Int = columns / 2 + 1
    val temp: Array[Array[Double]] = Array.ofDim[Double](n2d2, 2 * rows)
    val nthreads: Int = ConcurrencyUtils.numThreads
    if ((nthreads > 1) && useThreads && (rows >= nthreads) && (n2d2 - 2 >= nthreads)) {
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      var p: Int = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (i <- firstRow until lastRow) {
            fftColumns.realInverse2(a, i * columns, scale)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (r <- 0 until rows) {
        temp(0)(r) = a(r * columns)
      }
      fftRows.realInverseFull(temp(0), scale)
      p = (n2d2 - 2) / nthreads
      for (l <- 0 until nthreads) {
        val firstColumn: Int = 1 + l * p
        val lastColumn: Int = if (l == (nthreads - 1)) {
          n2d2 - 1
        }
        else {
          firstColumn + p
        }
        futures(l) = Future {
          for (c <- firstColumn until lastColumn) {
            val idx0: Int = 2 * c
            for (r <- 0 until rows) {
              val idx1: Int = 2 * r
              val idx2: Int = r * columns + idx0
              temp(c)(idx1) = a(idx2)
              temp(c)(idx1 + 1) = a(idx2 + 1)
            }
            fftRows.complexInverse(temp(c), scale)
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r * columns + 1)
        }
        fftRows.realInverseFull(temp(n2d2 - 1), scale)
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = r * columns
          val idx3: Int = n2d2 - 1
          temp(idx3)(idx1) = a(idx2 + 2 * idx3)
          temp(idx3)(idx1 + 1) = a(idx2 + 1)
        }
        fftRows.complexInverse(temp(n2d2 - 1), scale)
      }
      p = rows / nthreads
      for (l <- 0 until nthreads) {
        val firstRow: Int = l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx1: Int = 2 * r
            for (c <- 0 until n2d2) {
              val idx0: Int = 2 * c
              val idx2: Int = r * rowStride + idx0
              a(idx2) = temp(c)(idx1)
              a(idx2 + 1) = temp(c)(idx1 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
      for (l <- 0 until nthreads) {
        val firstRow: Int = 1 + l * p
        val lastRow: Int = if (l == (nthreads - 1)) {
          rows
        }
        else {
          firstRow + p
        }
        futures(l) = Future {
          for (r <- firstRow until lastRow) {
            val idx5: Int = r * rowStride
            val idx6: Int = (rows - r + 1) * rowStride
            for (c <- n2d2 until columns) {
              val idx1: Int = 2 * c
              val idx2: Int = 2 * (columns - c)
              a(idx1) = a(idx2)
              a(idx1 + 1) = -a(idx2 + 1)
              val idx3: Int = idx5 + idx1
              val idx4: Int = idx6 - idx1
              a(idx3) = a(idx4)
              a(idx3 + 1) = -a(idx4 + 1)
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (r <- 0 until rows) {
        fftColumns.realInverse2(a, r * columns, scale)
      }
      for (r <- 0 until rows) {
        temp(0)(r) = a(r * columns)
      }
      fftRows.realInverseFull(temp(0), scale)
      for (c <- 1 until n2d2 - 1) {
        val idx0: Int = 2 * c
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = r * columns + idx0
          temp(c)(idx1) = a(idx2)
          temp(c)(idx1 + 1) = a(idx2 + 1)
        }
        fftRows.complexInverse(temp(c), scale)
      }
      if ((columns % 2) == 0) {
        for (r <- 0 until rows) {
          temp(n2d2 - 1)(r) = a(r * columns + 1)
        }
        fftRows.realInverseFull(temp(n2d2 - 1), scale)
      }
      else {
        for (r <- 0 until rows) {
          val idx1: Int = 2 * r
          val idx2: Int = r * columns
          val idx3: Int = n2d2 - 1
          temp(idx3)(idx1) = a(idx2 + 2 * idx3)
          temp(idx3)(idx1 + 1) = a(idx2 + 1)
        }
        fftRows.complexInverse(temp(n2d2 - 1), scale)
      }
      for (r <- 0 until rows) {
        val idx1: Int = 2 * r
        for (c <- 0 until n2d2) {
          val idx0: Int = 2 * c
          val idx2: Int = r * rowStride + idx0
          a(idx2) = temp(c)(idx1)
          a(idx2 + 1) = temp(c)(idx1 + 1)
        }
      }
      for (r <- 1 until rows) {
        val idx5: Int = r * rowStride
        val idx6: Int = (rows - r + 1) * rowStride
        for (c <- n2d2 until columns) {
          val idx1: Int = 2 * c
          val idx2: Int = 2 * (columns - c)
          a(idx1) = a(idx2)
          a(idx1 + 1) = -a(idx2 + 1)
          val idx3: Int = idx5 + idx1
          val idx4: Int = idx6 - idx1
          a(idx3) = a(idx4)
          a(idx3 + 1) = -a(idx4 + 1)
        }
      }
    }
  }

  private def rdft2d_sub(isgn: Int, a: Array[Double]): Unit = {
    var n1h: Int = 0
    var j: Int = 0
    var xi: Double = .0
    var idx1: Int = 0
    var idx2: Int = 0
    n1h = rows >> 1
    if (isgn < 0) {
      for (i <- 1 until n1h) {
        j = rows - i
        idx1 = i * columns
        idx2 = j * columns
        xi = a(idx1) - a(idx2)
        a(idx1) += a(idx2)
        a(idx2) = xi
        xi = a(idx2 + 1) - a(idx1 + 1)
        a(idx1 + 1) += a(idx2 + 1)
        a(idx2 + 1) = xi
      }
    }
    else {
      for (i <- 1 until n1h) {
        j = rows - i
        idx1 = i * columns
        idx2 = j * columns
        a(idx2) = 0.5f * (a(idx1) - a(idx2))
        a(idx1) -= a(idx2)
        a(idx2 + 1) = 0.5f * (a(idx1 + 1) + a(idx2 + 1))
        a(idx1 + 1) -= a(idx2 + 1)
      }
    }
  }

  private def rdft2d_sub(isgn: Int, a: Array[Array[Double]]): Unit = {
    var n1h: Int = 0
    var j: Int = 0
    var xi: Double = .0
    n1h = rows >> 1
    if (isgn < 0) {
      for (i <- 1 until n1h) {
        j = rows - i
        xi = a(i)(0) - a(j)(0)
        a(i)(0) += a(j)(0)
        a(j)(0) = xi
        xi = a(j)(1) - a(i)(1)
        a(i)(1) += a(j)(1)
        a(j)(1) = xi
      }
    }
    else {
      for (i <- 1 until n1h) {
        j = rows - i
        a(j)(0) = 0.5f * (a(i)(0) - a(j)(0))
        a(i)(0) -= a(j)(0)
        a(j)(1) = 0.5f * (a(i)(1) + a(j)(1))
        a(i)(1) -= a(j)(1)
      }
    }
  }

  private def cdft2d_sub(isgn: Int, a: Array[Double], scale: Boolean): Unit = {
    var idx1: Int = 0
    var idx2: Int = 0
    var idx3: Int = 0
    var idx4: Int = 0
    var idx5: Int = 0
    var nt: Int = 8 * rows
    if (columns == 4) {
      nt >>= 1
    }
    else {
      if (columns < 4) {
        nt >>= 2
      }
    }
    val t: Array[Double] = new Array[Double](nt)
    if (isgn == -1) {
      if (columns > 4) {
        var c: Int = 0
        while ( {
          c < columns
        }) {
          for (r <- 0 until rows) {
            idx1 = r * columns + c
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            t(idx2) = a(idx1)
            t(idx2 + 1) = a(idx1 + 1)
            t(idx3) = a(idx1 + 2)
            t(idx3 + 1) = a(idx1 + 3)
            t(idx4) = a(idx1 + 4)
            t(idx4 + 1) = a(idx1 + 5)
            t(idx5) = a(idx1 + 6)
            t(idx5 + 1) = a(idx1 + 7)
          }
          fftRows.complexForward(t, 0)
          fftRows.complexForward(t, 2 * rows)
          fftRows.complexForward(t, 4 * rows)
          fftRows.complexForward(t, 6 * rows)
          for (r <- 0 until rows) {
            idx1 = r * columns + c
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            a(idx1) = t(idx2)
            a(idx1 + 1) = t(idx2 + 1)
            a(idx1 + 2) = t(idx3)
            a(idx1 + 3) = t(idx3 + 1)
            a(idx1 + 4) = t(idx4)
            a(idx1 + 5) = t(idx4 + 1)
            a(idx1 + 6) = t(idx5)
            a(idx1 + 7) = t(idx5 + 1)
          }

          c += 8
        }
      }
      else {
        if (columns == 4) {
          for (r <- 0 until rows) {
            idx1 = r * columns
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            t(idx2) = a(idx1)
            t(idx2 + 1) = a(idx1 + 1)
            t(idx3) = a(idx1 + 2)
            t(idx3 + 1) = a(idx1 + 3)
          }
          fftRows.complexForward(t, 0)
          fftRows.complexForward(t, 2 * rows)
          for (r <- 0 until rows) {
            idx1 = r * columns
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            a(idx1) = t(idx2)
            a(idx1 + 1) = t(idx2 + 1)
            a(idx1 + 2) = t(idx3)
            a(idx1 + 3) = t(idx3 + 1)
          }
        }
        else {
          if (columns == 2) {
            for (r <- 0 until rows) {
              idx1 = r * columns
              idx2 = 2 * r
              t(idx2) = a(idx1)
              t(idx2 + 1) = a(idx1 + 1)
            }
            fftRows.complexForward(t, 0)
            for (r <- 0 until rows) {
              idx1 = r * columns
              idx2 = 2 * r
              a(idx1) = t(idx2)
              a(idx1 + 1) = t(idx2 + 1)
            }
          }
        }
      }
    }
    else {
      if (columns > 4) {
        var c: Int = 0
        while ( {
          c < columns
        }) {
          for (r <- 0 until rows) {
            idx1 = r * columns + c
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            t(idx2) = a(idx1)
            t(idx2 + 1) = a(idx1 + 1)
            t(idx3) = a(idx1 + 2)
            t(idx3 + 1) = a(idx1 + 3)
            t(idx4) = a(idx1 + 4)
            t(idx4 + 1) = a(idx1 + 5)
            t(idx5) = a(idx1 + 6)
            t(idx5 + 1) = a(idx1 + 7)
          }
          fftRows.complexInverse(t, 0, scale)
          fftRows.complexInverse(t, 2 * rows, scale)
          fftRows.complexInverse(t, 4 * rows, scale)
          fftRows.complexInverse(t, 6 * rows, scale)
          for (r <- 0 until rows) {
            idx1 = r * columns + c
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            a(idx1) = t(idx2)
            a(idx1 + 1) = t(idx2 + 1)
            a(idx1 + 2) = t(idx3)
            a(idx1 + 3) = t(idx3 + 1)
            a(idx1 + 4) = t(idx4)
            a(idx1 + 5) = t(idx4 + 1)
            a(idx1 + 6) = t(idx5)
            a(idx1 + 7) = t(idx5 + 1)
          }

          c += 8
        }
      }
      else {
        if (columns == 4) {
          for (r <- 0 until rows) {
            idx1 = r * columns
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            t(idx2) = a(idx1)
            t(idx2 + 1) = a(idx1 + 1)
            t(idx3) = a(idx1 + 2)
            t(idx3 + 1) = a(idx1 + 3)
          }
          fftRows.complexInverse(t, 0, scale)
          fftRows.complexInverse(t, 2 * rows, scale)
          for (r <- 0 until rows) {
            idx1 = r * columns
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            a(idx1) = t(idx2)
            a(idx1 + 1) = t(idx2 + 1)
            a(idx1 + 2) = t(idx3)
            a(idx1 + 3) = t(idx3 + 1)
          }
        }
        else {
          if (columns == 2) {
            for (r <- 0 until rows) {
              idx1 = r * columns
              idx2 = 2 * r
              t(idx2) = a(idx1)
              t(idx2 + 1) = a(idx1 + 1)
            }
            fftRows.complexInverse(t, 0, scale)
            for (r <- 0 until rows) {
              idx1 = r * columns
              idx2 = 2 * r
              a(idx1) = t(idx2)
              a(idx1 + 1) = t(idx2 + 1)
            }
          }
        }
      }
    }
  }

  private def cdft2d_sub(isgn: Int, a: Array[Array[Double]], scale: Boolean): Unit = {
    var idx2: Int = 0
    var idx3: Int = 0
    var idx4: Int = 0
    var idx5: Int = 0
    var nt: Int = 8 * rows
    if (columns == 4) {
      nt >>= 1
    }
    else {
      if (columns < 4) {
        nt >>= 2
      }
    }
    val t: Array[Double] = new Array[Double](nt)
    if (isgn == -1) {
      if (columns > 4) {
        var c: Int = 0
        while ( {
          c < columns
        }) {
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            t(idx2) = a(r)(c)
            t(idx2 + 1) = a(r)(c + 1)
            t(idx3) = a(r)(c + 2)
            t(idx3 + 1) = a(r)(c + 3)
            t(idx4) = a(r)(c + 4)
            t(idx4 + 1) = a(r)(c + 5)
            t(idx5) = a(r)(c + 6)
            t(idx5 + 1) = a(r)(c + 7)
          }
          fftRows.complexForward(t, 0)
          fftRows.complexForward(t, 2 * rows)
          fftRows.complexForward(t, 4 * rows)
          fftRows.complexForward(t, 6 * rows)
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            a(r)(c) = t(idx2)
            a(r)(c + 1) = t(idx2 + 1)
            a(r)(c + 2) = t(idx3)
            a(r)(c + 3) = t(idx3 + 1)
            a(r)(c + 4) = t(idx4)
            a(r)(c + 5) = t(idx4 + 1)
            a(r)(c + 6) = t(idx5)
            a(r)(c + 7) = t(idx5 + 1)
          }

          c += 8
        }
      }
      else {
        if (columns == 4) {
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            t(idx2) = a(r)(0)
            t(idx2 + 1) = a(r)(1)
            t(idx3) = a(r)(2)
            t(idx3 + 1) = a(r)(3)
          }
          fftRows.complexForward(t, 0)
          fftRows.complexForward(t, 2 * rows)
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            a(r)(0) = t(idx2)
            a(r)(1) = t(idx2 + 1)
            a(r)(2) = t(idx3)
            a(r)(3) = t(idx3 + 1)
          }
        }
        else {
          if (columns == 2) {
            for (r <- 0 until rows) {
              idx2 = 2 * r
              t(idx2) = a(r)(0)
              t(idx2 + 1) = a(r)(1)
            }
            fftRows.complexForward(t, 0)
            for (r <- 0 until rows) {
              idx2 = 2 * r
              a(r)(0) = t(idx2)
              a(r)(1) = t(idx2 + 1)
            }
          }
        }
      }
    }
    else {
      if (columns > 4) {
        var c: Int = 0
        while ( {
          c < columns
        }) {
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            t(idx2) = a(r)(c)
            t(idx2 + 1) = a(r)(c + 1)
            t(idx3) = a(r)(c + 2)
            t(idx3 + 1) = a(r)(c + 3)
            t(idx4) = a(r)(c + 4)
            t(idx4 + 1) = a(r)(c + 5)
            t(idx5) = a(r)(c + 6)
            t(idx5 + 1) = a(r)(c + 7)
          }
          fftRows.complexInverse(t, 0, scale)
          fftRows.complexInverse(t, 2 * rows, scale)
          fftRows.complexInverse(t, 4 * rows, scale)
          fftRows.complexInverse(t, 6 * rows, scale)
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            idx4 = idx3 + 2 * rows
            idx5 = idx4 + 2 * rows
            a(r)(c) = t(idx2)
            a(r)(c + 1) = t(idx2 + 1)
            a(r)(c + 2) = t(idx3)
            a(r)(c + 3) = t(idx3 + 1)
            a(r)(c + 4) = t(idx4)
            a(r)(c + 5) = t(idx4 + 1)
            a(r)(c + 6) = t(idx5)
            a(r)(c + 7) = t(idx5 + 1)
          }

          c += 8
        }
      }
      else {
        if (columns == 4) {
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            t(idx2) = a(r)(0)
            t(idx2 + 1) = a(r)(1)
            t(idx3) = a(r)(2)
            t(idx3 + 1) = a(r)(3)
          }
          fftRows.complexInverse(t, 0, scale)
          fftRows.complexInverse(t, 2 * rows, scale)
          for (r <- 0 until rows) {
            idx2 = 2 * r
            idx3 = 2 * rows + 2 * r
            a(r)(0) = t(idx2)
            a(r)(1) = t(idx2 + 1)
            a(r)(2) = t(idx3)
            a(r)(3) = t(idx3 + 1)
          }
        }
        else {
          if (columns == 2) {
            for (r <- 0 until rows) {
              idx2 = 2 * r
              t(idx2) = a(r)(0)
              t(idx2 + 1) = a(r)(1)
            }
            fftRows.complexInverse(t, 0, scale)
            for (r <- 0 until rows) {
              idx2 = 2 * r
              a(r)(0) = t(idx2)
              a(r)(1) = t(idx2 + 1)
            }
          }
        }
      }
    }
  }

  private def xdft2d0_subth1(icr: Int, isgn: Int, a: Array[Double], scale: Boolean): Unit = {
    val nthreads: Int = min(rows, ConcurrencyUtils.numThreads)
    val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
    for (i <- 0 until nthreads) {
      val n0: Int = i
      futures(i) = Future {
        if (icr == 0) {
          if (isgn == -1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexForward(a, r * columns)

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexInverse(a, r * columns, scale)

              r += nthreads
            }
          }
        }
        else {
          if (isgn == 1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realForward(a, r * columns)

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realInverse(a, r * columns, scale)

              r += nthreads
            }
          }
        }
      }
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  private def xdft2d0_subth2(icr: Int, isgn: Int, a: Array[Double], scale: Boolean): Unit = {
    val nthreads: Int = min(rows, ConcurrencyUtils.numThreads)
    val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
    for (i <- 0 until nthreads) {
      val n0: Int = i
      futures(i) = Future {
        if (icr == 0) {
          if (isgn == -1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexForward(a, r * columns)

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexInverse(a, r * columns, scale)

              r += nthreads
            }
          }
        }
        else {
          if (isgn == 1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realForward(a, r * columns)

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realInverse2(a, r * columns, scale)

              r += nthreads
            }
          }
        }
      }
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  private def xdft2d0_subth1(icr: Int, isgn: Int, a: Array[Array[Double]], scale: Boolean): Unit = {
    val nthreads: Int = min(rows, ConcurrencyUtils.numThreads)
    val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
    for (i <- 0 until nthreads) {
      val n0: Int = i
      futures(i) = Future {
        if (icr == 0) {
          if (isgn == -1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexForward(a(r))

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexInverse(a(r), scale)

              r += nthreads
            }
          }
        }
        else {
          if (isgn == 1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realForward(a(r))

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realInverse(a(r), scale)

              r += nthreads
            }
          }
        }
      }
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  private def xdft2d0_subth2(icr: Int, isgn: Int, a: Array[Array[Double]], scale: Boolean): Unit = {
    val nthreads: Int = min(rows, ConcurrencyUtils.numThreads)
    val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
    for (i <- 0 until nthreads) {
      val n0: Int = i
      futures(i) = Future {
        if (icr == 0) {
          if (isgn == -1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexForward(a(r))

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.complexInverse(a(r), scale)

              r += nthreads
            }
          }
        }
        else {
          if (isgn == 1) {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realForward(a(r))

              r += nthreads
            }
          }
          else {
            var r: Int = n0
            while ( {
              r < rows
            }) {
              fftColumns.realInverse2(a(r), 0, scale)

              r += nthreads
            }
          }
        }
      }
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  private def cdft2d_subth(isgn: Int, a: Array[Double], scale: Boolean): Unit = {
    val nthread: Int = min(columns / 2, ConcurrencyUtils.numThreads)
    var nt: Int = 8 * rows
    if (columns == 4) {
      nt >>= 1
    }
    else {
      if (columns < 4) {
        nt >>= 2
      }
    }
    val ntf: Int = nt
    val futures: Array[Future[_]] = new Array[Future[_]](nthread)
    val nthreads: Int = nthread
    for (i <- 0 until nthread) {
      val n0: Int = i
      futures(i) = Future {
        var idx1: Int = 0
        var idx2: Int = 0
        var idx3: Int = 0
        var idx4: Int = 0
        var idx5: Int = 0
        val t: Array[Double] = new Array[Double](ntf)
        if (isgn == -1) {
          if (columns > 4 * nthreads) {
            var c: Int = 8 * n0
            while ( {
              c < columns
            }) {
              for (r <- 0 until rows) {
                idx1 = r * columns + c
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                t(idx2) = a(idx1)
                t(idx2 + 1) = a(idx1 + 1)
                t(idx3) = a(idx1 + 2)
                t(idx3 + 1) = a(idx1 + 3)
                t(idx4) = a(idx1 + 4)
                t(idx4 + 1) = a(idx1 + 5)
                t(idx5) = a(idx1 + 6)
                t(idx5 + 1) = a(idx1 + 7)
              }
              fftRows.complexForward(t, 0)
              fftRows.complexForward(t, 2 * rows)
              fftRows.complexForward(t, 4 * rows)
              fftRows.complexForward(t, 6 * rows)
              for (r <- 0 until rows) {
                idx1 = r * columns + c
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                a(idx1) = t(idx2)
                a(idx1 + 1) = t(idx2 + 1)
                a(idx1 + 2) = t(idx3)
                a(idx1 + 3) = t(idx3 + 1)
                a(idx1 + 4) = t(idx4)
                a(idx1 + 5) = t(idx4 + 1)
                a(idx1 + 6) = t(idx5)
                a(idx1 + 7) = t(idx5 + 1)
              }

              c += 8 * nthreads
            }
          }
          else {
            if (columns == 4 * nthreads) {
              for (r <- 0 until rows) {
                idx1 = r * columns + 4 * n0
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                t(idx2) = a(idx1)
                t(idx2 + 1) = a(idx1 + 1)
                t(idx3) = a(idx1 + 2)
                t(idx3 + 1) = a(idx1 + 3)
              }
              fftRows.complexForward(t, 0)
              fftRows.complexForward(t, 2 * rows)
              for (r <- 0 until rows) {
                idx1 = r * columns + 4 * n0
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                a(idx1) = t(idx2)
                a(idx1 + 1) = t(idx2 + 1)
                a(idx1 + 2) = t(idx3)
                a(idx1 + 3) = t(idx3 + 1)
              }
            }
            else {
              if (columns == 2 * nthreads) {
                for (r <- 0 until rows) {
                  idx1 = r * columns + 2 * n0
                  idx2 = 2 * r
                  t(idx2) = a(idx1)
                  t(idx2 + 1) = a(idx1 + 1)
                }
                fftRows.complexForward(t, 0)
                for (r <- 0 until rows) {
                  idx1 = r * columns + 2 * n0
                  idx2 = 2 * r
                  a(idx1) = t(idx2)
                  a(idx1 + 1) = t(idx2 + 1)
                }
              }
            }
          }
        }
        else {
          if (columns > 4 * nthreads) {
            var c: Int = 8 * n0
            while ( {
              c < columns
            }) {
              for (r <- 0 until rows) {
                idx1 = r * columns + c
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                t(idx2) = a(idx1)
                t(idx2 + 1) = a(idx1 + 1)
                t(idx3) = a(idx1 + 2)
                t(idx3 + 1) = a(idx1 + 3)
                t(idx4) = a(idx1 + 4)
                t(idx4 + 1) = a(idx1 + 5)
                t(idx5) = a(idx1 + 6)
                t(idx5 + 1) = a(idx1 + 7)
              }
              fftRows.complexInverse(t, 0, scale)
              fftRows.complexInverse(t, 2 * rows, scale)
              fftRows.complexInverse(t, 4 * rows, scale)
              fftRows.complexInverse(t, 6 * rows, scale)
              for (r <- 0 until rows) {
                idx1 = r * columns + c
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                a(idx1) = t(idx2)
                a(idx1 + 1) = t(idx2 + 1)
                a(idx1 + 2) = t(idx3)
                a(idx1 + 3) = t(idx3 + 1)
                a(idx1 + 4) = t(idx4)
                a(idx1 + 5) = t(idx4 + 1)
                a(idx1 + 6) = t(idx5)
                a(idx1 + 7) = t(idx5 + 1)
              }

              c += 8 * nthreads
            }
          }
          else {
            if (columns == 4 * nthreads) {
              for (r <- 0 until rows) {
                idx1 = r * columns + 4 * n0
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                t(idx2) = a(idx1)
                t(idx2 + 1) = a(idx1 + 1)
                t(idx3) = a(idx1 + 2)
                t(idx3 + 1) = a(idx1 + 3)
              }
              fftRows.complexInverse(t, 0, scale)
              fftRows.complexInverse(t, 2 * rows, scale)
              for (r <- 0 until rows) {
                idx1 = r * columns + 4 * n0
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                a(idx1) = t(idx2)
                a(idx1 + 1) = t(idx2 + 1)
                a(idx1 + 2) = t(idx3)
                a(idx1 + 3) = t(idx3 + 1)
              }
            }
            else {
              if (columns == 2 * nthreads) {
                for (r <- 0 until rows) {
                  idx1 = r * columns + 2 * n0
                  idx2 = 2 * r
                  t(idx2) = a(idx1)
                  t(idx2 + 1) = a(idx1 + 1)
                }
                fftRows.complexInverse(t, 0, scale)
                for (r <- 0 until rows) {
                  idx1 = r * columns + 2 * n0
                  idx2 = 2 * r
                  a(idx1) = t(idx2)
                  a(idx1 + 1) = t(idx2 + 1)
                }
              }
            }
          }
        }
      }
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  private def cdft2d_subth(isgn: Int, a: Array[Array[Double]], scale: Boolean): Unit = {
    val nthread: Int = min(columns / 2, ConcurrencyUtils.numThreads)
    var nt: Int = 8 * rows
    if (columns == 4) {
      nt >>= 1
    }
    else {
      if (columns < 4) {
        nt >>= 2
      }
    }
    val ntf: Int = nt
    val futures: Array[Future[_]] = new Array[Future[_]](nthread)
    val nthreads: Int = nthread
    for (i <- 0 until nthreads) {
      val n0: Int = i
      futures(i) = Future {
        var idx2: Int = 0
        var idx3: Int = 0
        var idx4: Int = 0
        var idx5: Int = 0
        val t: Array[Double] = new Array[Double](ntf)
        if (isgn == -1) {
          if (columns > 4 * nthreads) {
            var c: Int = 8 * n0
            while ( {
              c < columns
            }) {
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                t(idx2) = a(r)(c)
                t(idx2 + 1) = a(r)(c + 1)
                t(idx3) = a(r)(c + 2)
                t(idx3 + 1) = a(r)(c + 3)
                t(idx4) = a(r)(c + 4)
                t(idx4 + 1) = a(r)(c + 5)
                t(idx5) = a(r)(c + 6)
                t(idx5 + 1) = a(r)(c + 7)
              }
              fftRows.complexForward(t, 0)
              fftRows.complexForward(t, 2 * rows)
              fftRows.complexForward(t, 4 * rows)
              fftRows.complexForward(t, 6 * rows)
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                a(r)(c) = t(idx2)
                a(r)(c + 1) = t(idx2 + 1)
                a(r)(c + 2) = t(idx3)
                a(r)(c + 3) = t(idx3 + 1)
                a(r)(c + 4) = t(idx4)
                a(r)(c + 5) = t(idx4 + 1)
                a(r)(c + 6) = t(idx5)
                a(r)(c + 7) = t(idx5 + 1)
              }

              c += 8 * nthreads
            }
          }
          else {
            if (columns == 4 * nthreads) {
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                t(idx2) = a(r)(4 * n0)
                t(idx2 + 1) = a(r)(4 * n0 + 1)
                t(idx3) = a(r)(4 * n0 + 2)
                t(idx3 + 1) = a(r)(4 * n0 + 3)
              }
              fftRows.complexForward(t, 0)
              fftRows.complexForward(t, 2 * rows)
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                a(r)(4 * n0) = t(idx2)
                a(r)(4 * n0 + 1) = t(idx2 + 1)
                a(r)(4 * n0 + 2) = t(idx3)
                a(r)(4 * n0 + 3) = t(idx3 + 1)
              }
            }
            else {
              if (columns == 2 * nthreads) {
                for (r <- 0 until rows) {
                  idx2 = 2 * r
                  t(idx2) = a(r)(2 * n0)
                  t(idx2 + 1) = a(r)(2 * n0 + 1)
                }
                fftRows.complexForward(t, 0)
                for (r <- 0 until rows) {
                  idx2 = 2 * r
                  a(r)(2 * n0) = t(idx2)
                  a(r)(2 * n0 + 1) = t(idx2 + 1)
                }
              }
            }
          }
        }
        else {
          if (columns > 4 * nthreads) {
            var c: Int = 8 * n0
            while ( {
              c < columns
            }) {
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                t(idx2) = a(r)(c)
                t(idx2 + 1) = a(r)(c + 1)
                t(idx3) = a(r)(c + 2)
                t(idx3 + 1) = a(r)(c + 3)
                t(idx4) = a(r)(c + 4)
                t(idx4 + 1) = a(r)(c + 5)
                t(idx5) = a(r)(c + 6)
                t(idx5 + 1) = a(r)(c + 7)
              }
              fftRows.complexInverse(t, 0, scale)
              fftRows.complexInverse(t, 2 * rows, scale)
              fftRows.complexInverse(t, 4 * rows, scale)
              fftRows.complexInverse(t, 6 * rows, scale)
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                idx4 = idx3 + 2 * rows
                idx5 = idx4 + 2 * rows
                a(r)(c) = t(idx2)
                a(r)(c + 1) = t(idx2 + 1)
                a(r)(c + 2) = t(idx3)
                a(r)(c + 3) = t(idx3 + 1)
                a(r)(c + 4) = t(idx4)
                a(r)(c + 5) = t(idx4 + 1)
                a(r)(c + 6) = t(idx5)
                a(r)(c + 7) = t(idx5 + 1)
              }

              c += 8 * nthreads
            }
          }
          else {
            if (columns == 4 * nthreads) {
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                t(idx2) = a(r)(4 * n0)
                t(idx2 + 1) = a(r)(4 * n0 + 1)
                t(idx3) = a(r)(4 * n0 + 2)
                t(idx3 + 1) = a(r)(4 * n0 + 3)
              }
              fftRows.complexInverse(t, 0, scale)
              fftRows.complexInverse(t, 2 * rows, scale)
              for (r <- 0 until rows) {
                idx2 = 2 * r
                idx3 = 2 * rows + 2 * r
                a(r)(4 * n0) = t(idx2)
                a(r)(4 * n0 + 1) = t(idx2 + 1)
                a(r)(4 * n0 + 2) = t(idx3)
                a(r)(4 * n0 + 3) = t(idx3 + 1)
              }
            }
            else {
              if (columns == 2 * nthreads) {
                for (r <- 0 until rows) {
                  idx2 = 2 * r
                  t(idx2) = a(r)(2 * n0)
                  t(idx2 + 1) = a(r)(2 * n0 + 1)
                }
                fftRows.complexInverse(t, 0, scale)
                for (r <- 0 until rows) {
                  idx2 = 2 * r
                  a(r)(2 * n0) = t(idx2)
                  a(r)(2 * n0 + 1) = t(idx2 + 1)
                }
              }
            }
          }
        }
      }
    }
    ConcurrencyUtils.waitForCompletion(futures)
  }

  private def fillSymmetric(a: Array[Double]): Unit = {
    val twon2: Int = 2 * columns
    var idx1: Int = 0
    var idx2: Int = 0
    var idx3: Int = 0
    var idx4: Int = 0
    val n1d2: Int = rows / 2
    for (r <- rows - 1 to 1 by -1) {
      idx1 = r * columns
      idx2 = 2 * idx1
      var c: Int = 0
      while ( {
        c < columns
      }) {
        a(idx2 + c) = a(idx1 + c)
        a(idx1 + c) = 0
        a(idx2 + c + 1) = a(idx1 + c + 1)
        a(idx1 + c + 1) = 0

        c += 2
      }
    }
    val nthreads: Int = ConcurrencyUtils.numThreads
    if ((nthreads > 1) && useThreads && (n1d2 >= nthreads)) {
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      val l1k: Int = n1d2 / nthreads
      val newn2: Int = 2 * columns
      for (i <- 0 until nthreads) {
        val l1offa = if (i == 0) {
          i * l1k + 1
        } else {
          i * l1k
        }
        val l1stopa = i * l1k + l1k
        val l2offa  = i * l1k
        val l2stopa = if (i == nthreads - 1) {
          i * l1k + l1k + 1
        } else {
          i * l1k + l1k
        }
        futures(i) = Future {
          var idx1: Int = 0
          var idx2: Int = 0
          var idx3: Int = 0
          var idx4: Int = 0
          for (r <- l1offa until l1stopa) {
            idx1 = r * newn2
            idx2 = (rows - r) * newn2
            idx3 = idx1 + columns
            a(idx3) = a(idx2 + 1)
            a(idx3 + 1) = -a(idx2)
          }
          for (r <- l1offa until l1stopa) {
            idx1 = r * newn2
            idx3 = (rows - r + 1) * newn2
            var c: Int = columns + 2
            while ( {
              c < newn2
            }) {
              idx2 = idx3 - c
              idx4 = idx1 + c
              a(idx4) = a(idx2)
              a(idx4 + 1) = -a(idx2 + 1)

              c += 2
            }
          }
          for (r <- l2offa until l2stopa) {
            idx3 = ((rows - r) % rows) * newn2
            idx4 = r * newn2
            var c: Int = 0
            while ( {
              c < newn2
            }) {
              idx1 = idx3 + (newn2 - c) % newn2
              idx2 = idx4 + c
              a(idx1) = a(idx2)
              a(idx1 + 1) = -a(idx2 + 1)

              c += 2
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (r <- 1 until n1d2) {
        idx2 = r * twon2
        idx3 = (rows - r) * twon2
        a(idx2 + columns) = a(idx3 + 1)
        a(idx2 + columns + 1) = -a(idx3)
      }
      for (r <- 1 until n1d2) {
        idx2 = r * twon2
        idx3 = (rows - r + 1) * twon2
        var c: Int = columns + 2
        while ( {
          c < twon2
        }) {
          a(idx2 + c) = a(idx3 - c)
          a(idx2 + c + 1) = -a(idx3 - c + 1)

          c += 2
        }
      }
      for (r <- 0 to rows / 2) {
        idx1 = r * twon2
        idx4 = ((rows - r) % rows) * twon2
        var c: Int = 0
        while ( {
          c < twon2
        }) {
          idx2 = idx1 + c
          idx3 = idx4 + (twon2 - c) % twon2
          a(idx3) = a(idx2)
          a(idx3 + 1) = -a(idx2 + 1)

          c += 2
        }
      }
    }
    a(columns) = -a(1)
    a(1) = 0
    idx1 = n1d2 * twon2
    a(idx1 + columns) = -a(idx1 + 1)
    a(idx1 + 1) = 0
    a(idx1 + columns + 1) = 0
  }

  private def fillSymmetric(a: Array[Array[Double]]): Unit = {
    val newn2: Int = 2 * columns
    val n1d2: Int = rows / 2
    val nthreads: Int = ConcurrencyUtils.numThreads
    if ((nthreads > 1) && useThreads && (n1d2 >= nthreads)) {
      val futures: Array[Future[_]] = new Array[Future[_]](nthreads)
      val l1k: Int = n1d2 / nthreads
      for (i <- 0 until nthreads) {
        val l1offa = if (i == 0) {
           i * l1k + 1
        } else {
          i * l1k
        }
        val l1stopa = i * l1k + l1k
        val l2offa = i * l1k
        val l2stopa = if (i == nthreads - 1) {
          i * l1k + l1k + 1
        } else {
          i * l1k + l1k
        }
        futures(i) = Future {
          var idx1: Int = 0
          var idx2: Int = 0
          for (r <- l1offa until l1stopa) {
            idx1 = rows - r
            a(r)(columns) = a(idx1)(1)
            a(r)(columns + 1) = -a(idx1)(0)
          }
          for (r <- l1offa until l1stopa) {
            idx1 = rows - r
            var c: Int = columns + 2
            while ( {
              c < newn2
            }) {
              idx2 = newn2 - c
              a(r)(c) = a(idx1)(idx2)
              a(r)(c + 1) = -a(idx1)(idx2 + 1)

              c += 2
            }
          }
          for (r <- l2offa until l2stopa) {
            idx1 = (rows - r) % rows
            var c: Int = 0
            while ( {
              c < newn2
            }) {
              idx2 = (newn2 - c) % newn2
              a(idx1)(idx2) = a(r)(c)
              a(idx1)(idx2 + 1) = -a(r)(c + 1)

              c = c + 2
            }
          }
        }
      }
      ConcurrencyUtils.waitForCompletion(futures)
    }
    else {
      for (r <- 1 until n1d2) {
        val idx1: Int = rows - r
        a(r)(columns) = a(idx1)(1)
        a(r)(columns + 1) = -a(idx1)(0)
      }
      for (r <- 1 until n1d2) {
        val idx1: Int = rows - r
        var c: Int = columns + 2
        while ( {
          c < newn2
        }) {
          val idx2: Int = newn2 - c
          a(r)(c) = a(idx1)(idx2)
          a(r)(c + 1) = -a(idx1)(idx2 + 1)

          c += 2
        }
      }
      for (r <- 0 to rows / 2) {
        val idx1: Int = (rows - r) % rows
        var c: Int = 0
        while ( {
          c < newn2
        }) {
          val idx2: Int = (newn2 - c) % newn2
          a(idx1)(idx2) = a(r)(c)
          a(idx1)(idx2 + 1) = -a(r)(c + 1)

          c += 2
        }
      }
    }
    a(0)(columns) = -a(0)(1)
    a(0)(1) = 0
    a(n1d2)(columns) = -a(n1d2)(1)
    a(n1d2)(1) = 0
    a(n1d2)(columns + 1) = 0
  }

}
