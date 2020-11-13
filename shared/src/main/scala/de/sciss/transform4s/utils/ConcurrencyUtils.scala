/* ***** BEGIN LICENSE BLOCK *****
 * JLargeArrays
 * Copyright (C) 2013 onward University of Warsaw, ICM
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

import java.lang.Math.max

import scala.concurrent.ExecutionContext

/**
 * Concurrency utilities.
 * <p>
 *
 * @author Piotr Wendykier (p.wendykier@icm.edu.pl)
 */
object ConcurrencyUtils extends ConcurrencyUtilsPlatform {
  private var nThreads: Int = numProcessors

  private var _concurrentThreshold: Int = 100000

  /**
   * Returns the minimum length of array for which multiple threads are used.
   * <p>
   *
   * @return the minimum length of array for which multiple threads are used
   */
  def concurrentThreshold: Int = _concurrentThreshold

  /**
   * Sets the minimum length of an array for which multiple threads are used.
   * <p>
   *
   * @param value minimum length of an array for which multiple threads are used
   */
  def concurrentThreshold_=(value: Int): Unit =
    _concurrentThreshold = max(1, value)

  /**
   * Returns the number of available processors.
   * <p>
   *
   * @return number of available processors
   */
  def numProcessors: Int = Runtime.getRuntime.availableProcessors

  /**
   * Returns the current number of threads.
   * <p>
   *
   * @return the current number of threads.
   */
  def numThreads: Int = nThreads

  /**
   * Sets the number of threads.
   * <p>
   *
   * @param n new value of threads
   */
  def numThreads_=(n: Int): Unit =
    nThreads = n

  /**
   * Sets the pool of threads.
   * <p>
   *
   * @param ec pool of threads
   */
  def executionContext_=(ec: ExecutionContext): Unit =
    _executionContext = ec

  /**
   * Returns the pool of threads.
   * <p>
   *
   * @return pool of threads
   */
  implicit def executionContext: ExecutionContext = _executionContext
}
