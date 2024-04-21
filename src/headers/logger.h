//
// Created by mafu on 2024/4/7.
//
#pragma once
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <thread>
#include <sstream>
template <>
struct fmt::formatter<std::thread::id> {
    static constexpr auto parse(const format_parse_context & ctx) -> decltype(ctx.begin()) {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const std::thread::id& tid, FormatContext& ctx) const -> decltype(ctx.out()) {
        std::ostringstream oss;
        oss << tid;
        return fmt::format_to(ctx.out(), "{}", oss.str());
    }
};
namespace Coal {
    class Logger {
    public:

        Logger(const Logger&) = delete;
        Logger& operator=(const Logger&) = delete;

        static Logger& getInstance() {
            static Logger instance;
            return instance;
        }

        void init(const std::string& logFilePath) {
            const auto console_sink =
                    std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            const auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logFilePath, true);

            console_sink->set_pattern("[%^%L%$] [%Y-%m-%d %H:%M:%S] [thread %t] %v");
            file_sink->set_pattern("[%^%L%$] [%Y-%m-%d %H:%M:%S] [thread %t] %v");

            std::vector<spdlog::sink_ptr> sinks {console_sink, file_sink};
            logger = std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());
            logger->set_level(spdlog::level::info);
            logger->flush_on(spdlog::level::info);
        }

        std::shared_ptr<spdlog::logger> get() {
            return logger;
        }

    private:
        std::shared_ptr<spdlog::logger> logger;

        Logger() = default;

    };

}// namespace Coal
