package de.sciss.transform4s.utils

import java.util.concurrent.{Executors, ThreadFactory}

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}

trait ConcurrencyUtilsPlatform {
  private val DEFAULT_THREAD_POOL = Executors.newCachedThreadPool(
    new CustomThreadFactory(new CustomExceptionHandler))

  protected final var _executionContext: ExecutionContext =
    ExecutionContext.fromExecutor(DEFAULT_THREAD_POOL)

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
}
