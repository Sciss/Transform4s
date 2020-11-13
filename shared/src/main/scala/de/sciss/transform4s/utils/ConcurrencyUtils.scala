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
import java.util.concurrent.{Executors, ThreadFactory}

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}

/**
 * Concurrency utilities.
 * <p>
 *
 * @author Piotr Wendykier (p.wendykier@icm.edu.pl)
 */
object ConcurrencyUtils {
  private val DEFAULT_THREAD_POOL = Executors.newCachedThreadPool(
    new ConcurrencyUtils.CustomThreadFactory(new ConcurrencyUtils.CustomExceptionHandler))

  private var _executionContext: ExecutionContext /*Executor*/ =
    ExecutionContext.fromExecutor(DEFAULT_THREAD_POOL)

  private var nThreads: Int = numProcessors

  private var _concurrentThreshold: Long = 100000L

  private class CustomExceptionHandler extends Thread.UncaughtExceptionHandler {
    override def uncaughtException(t: Thread, e: Throwable): Unit = {
      e.printStackTrace()
    }
  }

  private object CustomThreadFactory {
    private val DEFAULT_FACTORY = Executors.defaultThreadFactory
  }

  private class CustomThreadFactory(val handler: Thread.UncaughtExceptionHandler) extends ThreadFactory {
    override def newThread(r: Runnable): Thread = {
      val t = CustomThreadFactory.DEFAULT_FACTORY.newThread(r)
      t.setUncaughtExceptionHandler(handler)
      t
    }
  }

  /**
   * Returns the minimum length of array for which multiple threads are used.
   * <p>
   *
   * @return the minimum length of array for which multiple threads are used
   */
  def concurrentThreshold: Long = _concurrentThreshold

  /**
   * Sets the minimum length of an array for which multiple threads are used.
   * <p>
   *
   * @param value minimum length of an array for which multiple threads are used
   */
  def concurrentThreshold_=(value: Long): Unit =
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

//  /**
//   * Submits a value-returning task for execution and returns a Future
//   * representing the pending results of the task.
//   * <p>
//   *
//   * @tparam A  type
//   * @param task task for execution
//   *             <p>
//   * @return handle to the task submitted for execution
//   */
//  def submit[A](task: => A): Future[A] =
//    Future(task)(_executionContext)


//  /**
//   * Submits a Runnable task for execution and returns a Future representing that task.
//   * <p>
//   *
//   * @param task task for execution
//   *             <p>
//   * @return handle to the task submitted for execution
//   */
//  def submit(task: Runnable): Future[_] = {
//    ConcurrencyUtils.executionContext.execute(task)
//  }

  /**
   * Waits for all threads to complete computation.
   * <p>
   *
   * @param futures list of handles to the tasks
   *                <p>
   * @throws ExecutionException   if the computation threw an exception
   * @throws InterruptedException if the current thread was interrupted while waiting
   */
  def waitForCompletion(futures: Array[Future[_]]): Unit = {
    var i = 0
    while (i < futures.length) {
      Await.result(futures(i), Duration.Inf)
      i += 1
    }
  }

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

//  /**
//   * Shutdowns all submitted tasks.
//   */
//  def shutdownThreadPoolAndAwaitTermination(): Unit = {
//    ConcurrencyUtils._executionContext.shutdown() // Disable new tasks from being submitted
//
//    try // Wait a while for existing tasks to terminate
//      if (!ConcurrencyUtils._executionContext.awaitTermination(60, TimeUnit.SECONDS)) {
//        ConcurrencyUtils._executionContext.shutdownNow // Cancel currently executing tasks
//
//        // Wait a while for tasks to respond to being cancelled
//        if (!ConcurrencyUtils._executionContext.awaitTermination(60, TimeUnit.SECONDS)) System.err.println("Pool did not terminate")
//      }
//    catch {
//      case _: InterruptedException =>
//        // (Re-)Cancel if current thread also interrupted
//        ConcurrencyUtils._executionContext.shutdownNow
//        // Preserve interrupt status
//        Thread.currentThread.interrupt()
//    }
//  }
}
