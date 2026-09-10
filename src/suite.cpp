#include "quake/suite.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace quake {

SuiteRunResult run_record_suite_parallel(const NonlinearDynamicModel& model,
                                         const std::vector<std::vector<double>>& motions,
                                         double dt,int workers,double tolerance,int max_iterations){
    if(dt<=0.0||workers<=0||max_iterations<=0)throw std::invalid_argument("invalid suite options");
    SuiteRunResult out;out.records.resize(motions.size());if(motions.empty()){out.workers=0;return out;}
    workers=std::max(1,std::min(workers,static_cast<int>(motions.size())));out.workers=workers;out.setup_factorizations=static_cast<std::size_t>(workers);
    std::atomic<std::size_t> next{0};std::exception_ptr failure;std::mutex failure_mutex;
    const auto start=std::chrono::steady_clock::now();std::vector<std::thread> threads;threads.reserve(static_cast<std::size_t>(workers));
    for(int w=0;w<workers;++w)threads.emplace_back([&]{
        try{
            PreparedWoodburyNewmark prepared(model,dt);
            while(true){const std::size_t i=next.fetch_add(1);if(i>=motions.size())break;out.records[i]=prepared.run(motions[i],tolerance,max_iterations);}
        }catch(...){std::lock_guard<std::mutex> lock(failure_mutex);if(!failure)failure=std::current_exception();next.store(motions.size());}
    });
    for(auto& t:threads) t.join();
    out.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if(failure) std::rethrow_exception(failure);
    return out;
}

} // namespace quake
