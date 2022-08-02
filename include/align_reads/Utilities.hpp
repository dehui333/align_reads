#ifndef ALIGN_READS_UTILITIES_HPP_
#define ALIGN_READS_UTILITIES_HPP_

#include "thread_pool/thread_pool.hpp"

namespace align_reads
{
    template <typename R>
    class Futures
    {
    public:
        Futures<R>(std::shared_ptr<thread_pool::ThreadPool> &pool, size_t reserve = 0) : pool(pool)
        {
            futures.reserve(reserve);
        }

        template <typename F, typename... Ts>
        void add_inputs(F &&routine, Ts &&...params)
        {
            futures.emplace_back(pool->Submit(routine, params...));
        }
        std::vector<R> get()
        {
            std::vector<R> results;
            results.reserve(futures.size());
            for (auto &f : futures)
            {
                results.push_back(f.get());
            }
            futures.clear();
            return results;
        }

        void finish()
        {
            for (auto &f : futures)
            {
                f.get();
            }
            futures.clear();
        }

    private:
        std::vector<std::future<R>> futures;
        std::shared_ptr<thread_pool::ThreadPool> pool;
    };
}

#endif // ALIGN_READS_UTILITIES_HPP_