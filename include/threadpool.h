#pragma once
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

/* 一个简单的例子
#include "threadpool.h"
#include <iostream>
#include <chrono>
#include <vector>

// 假设的任务函数：模拟每个线程处理一个文件
void processFastaChunk(int id, int duration_ms) {
	std::this_thread::sleep_for(std::chrono::milliseconds(duration_ms));  // 模拟耗时任务
	std::cout << "任务 " << id << " 完成，耗时 " << duration_ms << "ms\n";
}

int main() {
	// 创建一个包含 4 个线程的线程池
	ThreadPool pool(4);

	// 存储 future 对象，用于获取返回值或等待完成
	std::vector<std::future<void>> futures;

	// 启动多个任务
	for (int i = 0; i < 10; ++i) {
		int duration = 100 + (i % 5) * 100; // 每个任务持续 100~500ms
		futures.emplace_back(
			pool.enqueue(processFastaChunk, i, duration)
		);
	}

	// 等待所有任务完成（可选，若不关心结果可以省略）
	for (auto& fut : futures) {
		fut.get();  // 等待并捕获异常
	}

	// 或者使用 ThreadPool 提供的等待方法
	// pool.waitAllTasksDone();

	std::cout << "所有任务完成。\n";
	return 0;
}

*/

// modified from https://github.com/progschj/ThreadPool.git
class ThreadPool {
public:
	ThreadPool(size_t);
	template<class F, class... Args>
	auto enqueue(F&& f, Args&&... args)
		-> std::future<typename std::result_of<F(Args...)>::type>;
	void waitAllTasksDone();
	~ThreadPool();
private:
	// need to keep track of threads so we can join them
	std::vector< std::thread > workers;
	// the task queue
	std::queue< std::function<void()> > tasks;

	// synchronization
	std::mutex queue_mutex;
	std::condition_variable condition;
	bool stop;

	std::condition_variable all_tasks_done;
	int tasks_count = 0;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
	: stop(false)
{
	for (size_t i = 0; i < threads; ++i)
		workers.emplace_back(
			[this]
			{
				for (;;)
				{
					std::function<void()> task;

					{
						std::unique_lock<std::mutex> lock(this->queue_mutex);
						this->condition.wait(lock,
							[this] { return this->stop || !this->tasks.empty(); });
						if (this->stop && this->tasks.empty())
							return;
						task = std::move(this->tasks.front());
						this->tasks.pop();
					}

					task();
				}
			}
		);
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
-> std::future<typename std::result_of<F(Args...)>::type>
{
	using return_type = typename std::result_of<F(Args...)>::type;

	auto task = std::make_shared< std::packaged_task<return_type()> >(
		std::bind(std::forward<F>(f), std::forward<Args>(args)...)
	);

	std::future<return_type> res = task->get_future();
	{
		std::unique_lock<std::mutex> lock(queue_mutex);

		// don't allow enqueueing after stopping the pool
		if (stop)
			throw std::runtime_error("enqueue on stopped ThreadPool");

		tasks_count++;

		tasks.emplace([task, this]() {
			(*task)();
			std::unique_lock<std::mutex> lock(this->queue_mutex);
			tasks_count--;
			if (tasks_count == 0) {
				all_tasks_done.notify_one();
			}
			});
	}
	condition.notify_one();
	return res;
}

inline void ThreadPool::waitAllTasksDone() {
	std::unique_lock<std::mutex> lock(queue_mutex);
	all_tasks_done.wait(lock, [this] { return tasks.empty() && tasks_count == 0; });
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
	{
		std::unique_lock<std::mutex> lock(queue_mutex);
		stop = true;
	}
	condition.notify_all();
	for (std::thread& worker : workers)
		worker.join();
}